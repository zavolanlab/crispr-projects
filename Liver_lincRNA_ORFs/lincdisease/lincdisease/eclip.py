"""
ENCODE HepG2 eCLIP features of the mechanism classes M4 (nuclear trans-acting RNA) and M5 (cytoplasmic
trans-acting RNA) (LOGIC.md, Step 4).

Peaks: reproducible (IDR, replicates 1+2) eCLIP peaks of every HepG2 RBP (GRCh38 narrowPeak; ENCODE file list in
ECLIP_FILES), kept at log2 fold enrichment over the size-matched input >= 3 and -log10 p >= 3 (Van Nostrand et al.
2020 thresholds). A feature is bound by an RBP if a peak on the feature's strand overlaps its exons minus the
exons of other GENCODE genes on the same strand (stranded_regions; eCLIP is stranded, so antisense genes are not
subtracted as they are for the unstranded public RNA-seq).

RBP groups, from the immunofluorescence localisation in HepG2 of Van Nostrand et al. 2020 (Nature 583:711,
Supplementary Data 1): nuclear = signal in a nuclear compartment (nuclei, nucleolus, speckles, PML or Cajal
bodies) and none in a cytoplasmic one (cytoplasm, mitochondria, Golgi, P bodies, ER); cytoplasmic = the reverse.
RBPs seen in both, or not imaged, belong to neither group (but count in the total).

Statistics per feature (enrichments missing if no RBP binds it; mechanism.features sets them to 0, no evidence,
for loci expressed in HepG2):
  eclip_n_rbps            number of RBPs bound (all groups)
  eclip_rbps_per_kb       eclip_n_rbps per kb of quantifiable sequence
  eclip_<group>_enrichment  log2 [ ((n_group + f_group) / (n_all + 1)) / f_group ], f_group = group size / number
                          of RBPs: over-representation of the group among the RBPs bound, shrunk toward 0 by one
                          pseudo-RBP split by the group fractions; 0 without binding, and the further from 0 the
                          more RBPs support the over- or under-representation.
"""
import gzip
import os

import numpy as np
import pandas as pd

from . import quant

ENCODE_SEARCH = ("https://www.encodeproject.org/search/?type=File&assay_title=eCLIP&biosample_ontology.term_name="
                 "HepG2&file_format=bed&assembly=GRCh38&status=released&format=json&limit=all"
                 "&field=accession&field=output_type&field=biological_replicates&field=target.label&field=dataset")
ENCODE_BED = "https://www.encodeproject.org/files/{0}/@@download/{0}.bed.gz"
RBP_ANNOTATION = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2077-3/MediaObjects/"
                  "41586_2020_2077_MOESM3_ESM.xlsx")
NUCLEAR = ["Nuclei", "Nucleolus", "Speckles", "PML bodies", "Cajal bodies"]
CYTOPLASMIC = ["Cytoplasm", "Mitochondria", "Golgi", "P bodies", "ER"]
MIN_LOG2FC = 3.0
MIN_MLOG10P = 3.0


def idr_files(search_json):
    """(accession, RBP) of the reproducible peak files (biological replicates 1 and 2) of an ENCODE search."""
    rows = [(f["accession"], f["target"]["label"]) for f in search_json["@graph"]
            if f.get("output_type") == "peaks" and sorted(f.get("biological_replicates", [])) == [1, 2]]
    return pd.DataFrame(rows, columns=["accession", "rbp"])


def rbp_groups(annotation_xlsx):
    """Series RBP -> 'nuclear' / 'cytoplasmic' / 'both' / '' (not imaged) from Van Nostrand 2020 Suppl. Data 1."""
    x = pd.read_excel(annotation_xlsx, header=None)
    cols = x.iloc[1].tolist()
    cols[0], cols[1] = "rbp", "gene_id"
    x = x.iloc[2:].set_axis(cols, axis=1).set_index("rbp")
    loc = x[NUCLEAR + CYTOPLASMIC].apply(pd.to_numeric, errors="coerce")
    return localisation_group(loc)


def localisation_group(loc):
    imaged = loc.notna().any(axis=1)
    nuc = loc[NUCLEAR].fillna(0).sum(axis=1) > 0
    cyt = loc[CYTOPLASMIC].fillna(0).sum(axis=1) > 0
    g = np.select([nuc & ~cyt, cyt & ~nuc, nuc & cyt], ["nuclear", "cytoplasmic", "both"], "")
    return pd.Series(np.where(imaged, g, ""), index=loc.index)


def read_peaks(path, min_log2fc=MIN_LOG2FC, min_mlog10p=MIN_MLOG10P):
    """chrom/start0/end1/strand of the peaks passing the enrichment thresholds."""
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 2, 5, 6, 7],
                         names=["chrom", "start0", "end1", "strand", "log2fc", "mlog10p"])
    df = df[(df["log2fc"] >= min_log2fc) & (df["mlog10p"] >= min_mlog10p)]
    return df[["chrom", "start0", "end1", "strand"]]


def stranded_regions(lincrna_exons, gencode_exons, loci, control_genes):
    """
    Exonic regions for stranded data: the union of the feature's exons minus exons of other GENCODE genes on the
    SAME strand (antisense genes do not compete for stranded reads; quant.build_regions removes both strands
    for unstranded RNA-seq, which e.g. leaves 11% of MALAT1). Returns (regions, strands).
    """
    lx = lincrna_exons.rename(columns={"gene_id": "feature_id"})
    lx = lx[lx["feature_id"].isin(loci.index)]
    ctrl = gencode_exons[gencode_exons["gene_name"].isin(control_genes)]
    cx = ctrl.assign(feature_id="ctrl:" + ctrl["gene_name"])
    ex = pd.concat([lx, cx], ignore_index=True)[["feature_id", "chrom", "start0", "end1", "strand"]]
    strands = ex.groupby("feature_id")["strand"].first()
    exons = quant.merge_intervals(ex.drop(columns="strand"))
    own = {lid: {i.split(".")[0] for i in quant._split_ids(r.get("name_source_gene_ids"))
                 + quant._split_ids(r.get("lncRNA_gene_ids"))} for lid, r in loci.iterrows()}
    for name in control_genes:
        own["ctrl:" + name] = set(ctrl.loc[ctrl["gene_name"] == name, "gene_id"].str.split(".").str[0])
    g = gencode_exons[["chrom", "start0", "end1", "strand"]].assign(
        gene_id=gencode_exons["gene_id"].str.split(".").str[0]).drop_duplicates()
    g_by = dict(tuple(g.groupby(["chrom", "strand"])))
    out = []
    for (fid, chrom), sub in exons.groupby(["feature_id", "chrom"]):
        b = g_by.get((chrom, strands[fid]), g.iloc[:0])
        b = b[(b["start0"] < sub["end1"].max()) & (b["end1"] > sub["start0"].min())]
        out.append(quant.subtract_intervals(sub, b[~b["gene_id"].isin(own.get(fid, set()))]))
    return quant.merge_intervals(pd.concat(out, ignore_index=True)), strands


def bound_features(regions, strands, peaks):
    """Set of feature ids with a peak on their strand overlapping one of their regions."""
    reg = regions.assign(strand=regions["feature_id"].map(strands)).dropna(subset=["strand"])
    hit = set()
    for (chrom, strand), r in reg.groupby(["chrom", "strand"]):
        p = peaks[(peaks["chrom"] == chrom) & (peaks["strand"] == strand)].sort_values("start0")
        if p.empty:
            continue
        ps, pe = p["start0"].to_numpy(), p["end1"].to_numpy()
        longest = int((pe - ps).max())
        # candidate peaks start in [region start - longest peak, region end)
        lo = np.searchsorted(ps, r["start0"].to_numpy() - longest, "left")
        hi = np.searchsorted(ps, r["end1"].to_numpy(), "left")
        starts = r["start0"].to_numpy()
        ok = [bool(np.any(pe[a:b] > s)) for a, b, s in zip(lo, hi, starts)]
        hit.update(r["feature_id"].to_numpy()[ok])
    return hit


def features(regions, strands, peak_sets, groups):
    """
    `peak_sets`: dict RBP -> peaks (read_peaks); `groups`: Series RBP -> group. One row per feature in `regions`.
    """
    fids = pd.Index(regions["feature_id"].unique())
    bound = pd.DataFrame(False, index=fids, columns=sorted(peak_sets))
    for rbp, peaks in peak_sets.items():
        b = bound_features(regions, strands, peaks)
        bound.loc[bound.index.isin(b), rbp] = True
    g = groups.reindex(bound.columns).fillna("")
    kb = (regions["end1"] - regions["start0"]).groupby(regions["feature_id"]).sum().reindex(fids) / 1000
    out = pd.DataFrame(index=fids)
    n_all = bound.sum(axis=1)
    out["eclip_n_rbps"] = n_all
    out["eclip_rbps_per_kb"] = n_all / kb
    for grp in ("nuclear", "cytoplasmic"):
        cols = g.index[g == grp]
        f = len(cols) / len(g)
        n = bound[cols].sum(axis=1)
        out[f"eclip_n_{grp}"] = n
        out[f"eclip_{grp}_enrichment"] = np.log2(((n + f) / (n_all + 1)) / f)
        out[f"eclip_{grp}_rbps"] = bound[cols].apply(lambda r: ",".join(r.index[r]), axis=1)
    out.loc[n_all == 0, ["eclip_rbps_per_kb", "eclip_nuclear_enrichment", "eclip_cytoplasmic_enrichment"]] = np.nan
    out["eclip_n_rbps_tested"] = len(g)
    return out


def load_peak_sets(files, data_dir):
    """files: DataFrame accession/rbp (idr_files). Replicated RBPs (none in HepG2 so far) are pooled."""
    sets = {}
    for acc, rbp in files[["accession", "rbp"]].itertuples(index=False):
        p = read_peaks(os.path.join(data_dir, f"{acc}.bed.gz"))
        sets[rbp] = pd.concat([sets[rbp], p]) if rbp in sets else p
    return sets
