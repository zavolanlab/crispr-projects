"""
Quantification of lincRNA loci and control genes in patient and tissue RNA-seq (LOGIC.md, Step 2) from
recount3 base-level coverage bigWigs.

Regions. Public RNA-seq (TCGA, GTEx, most SRA cohorts) is quantified unstranded, so a locus is measured on
its quantifiable exonic positions: the union of the exons of its isoforms minus the exons of every GENCODE
gene on either strand that is not the locus itself (the locus' own GENCODE lncRNA genes are kept). Control
genes are measured the same way on the union of their GENCODE exons. Loci with fewer than
`min_quant_bases` quantifiable positions are flagged as not quantifiable (QC flag `tcga_quantifiable`).

Values. For each sample, the summed base coverage over the region is divided by the number of positions and
by the sample's total base coverage (recount3 `bc_auc.all_reads_all_bases`) and multiplied by 1e9
(`ncov`, an RPKM-like unit independent of read length).

The bigWigs are streamed: downloaded one at a time into a scratch directory, summed over the regions and
deleted. Sums are appended to <out>/sums.<worker>.tsv (one row per sample), so an interrupted run resumes
where it stopped.
"""
import gzip
import os
import re
import shutil
import tempfile
import urllib.request

import numpy as np
import pandas as pd

RECOUNT3 = "https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/"


# ---------------------------------------------------------------- intervals

def merge_intervals(df, by=("feature_id", "chrom")):
    """Merge overlapping/adjacent intervals within each group; returns columns by + start0/end1."""
    out = []
    for key, g in df.sort_values(list(by) + ["start0"]).groupby(list(by), sort=False):
        s = g["start0"].to_numpy()
        e = g["end1"].to_numpy()
        cur_s, cur_e = s[0], e[0]
        for a, b in zip(s[1:], e[1:]):
            if a <= cur_e:
                cur_e = max(cur_e, b)
            else:
                out.append((*key, cur_s, cur_e))
                cur_s, cur_e = a, b
        out.append((*key, cur_s, cur_e))
    return pd.DataFrame(out, columns=list(by) + ["start0", "end1"])


def subtract_intervals(a, b):
    """Positions of the intervals `a` (feature_id/chrom/start0/end1) not covered by `b` (chrom/start0/end1)."""
    out = []
    b_by_chrom = {c: g.sort_values("start0") for c, g in b.groupby("chrom")}
    for row in a.itertuples(index=False):
        g = b_by_chrom.get(row.chrom)
        pieces = [(row.start0, row.end1)]
        if g is not None:
            hit = g[(g["start0"] < row.end1) & (g["end1"] > row.start0)]
            for bs, be in zip(hit["start0"], hit["end1"]):
                nxt = []
                for s, e in pieces:
                    if be <= s or bs >= e:
                        nxt.append((s, e))
                        continue
                    if s < bs:
                        nxt.append((s, bs))
                    if be < e:
                        nxt.append((be, e))
                pieces = nxt
        out.extend((row.feature_id, row.chrom, s, e) for s, e in pieces if e > s)
    return pd.DataFrame(out, columns=["feature_id", "chrom", "start0", "end1"])


# ---------------------------------------------------------------- annotation

def _attr(series, key):
    return series.str.extract(f'{key} "([^"]+)"', expand=False)


def read_gtf_exons(path):
    """Exons of a GTF: chrom/start0/end1/strand/gene_id/gene_name/gene_type/transcript_id."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                         names=["chrom", "feature", "start", "end", "strand", "attr"], dtype={0: str})
    df = df[df["feature"] == "exon"]
    return pd.DataFrame({
        "chrom": df["chrom"], "start0": df["start"] - 1, "end1": df["end"], "strand": df["strand"],
        "gene_id": _attr(df["attr"], "gene_id"), "gene_name": _attr(df["attr"], "gene_name"),
        "gene_type": _attr(df["attr"], "gene_type"), "transcript_id": _attr(df["attr"], "transcript_id"),
    }).reset_index(drop=True)


def _split_ids(value):
    if not isinstance(value, str) or not value.strip():
        return []
    return [v.strip() for v in value.split(",") if v.strip()]


def build_regions(lincrna_exons, gencode_exons, loci, control_genes, min_quant_bases=100):
    """
    Quantifiable regions of every lincRNA locus (feature_id = locus_id) and control gene
    (feature_id = 'ctrl:<gene_name>'). `loci` is inputs.loci() (index locus_id) and gives the GENCODE ids of
    the locus itself (name_source_gene_ids, lncRNA_gene_ids). Returns (regions, features).
    """
    lx = lincrna_exons.rename(columns={"gene_id": "feature_id"})[["feature_id", "chrom", "start0", "end1"]]
    lx = lx[lx["feature_id"].isin(loci.index)]
    ctrl = gencode_exons[gencode_exons["gene_name"].isin(control_genes)]
    missing = set(control_genes) - set(ctrl["gene_name"])
    if missing:
        raise ValueError(f"control genes not in the annotation: {sorted(missing)}")
    cx = ctrl.assign(feature_id="ctrl:" + ctrl["gene_name"])[["feature_id", "chrom", "start0", "end1"]]
    exons = merge_intervals(pd.concat([lx, cx], ignore_index=True))

    # GENCODE exons of the feature's own genes are not subtracted from it
    own = {}
    for locus_id, row in loci.iterrows():
        ids = _split_ids(row.get("name_source_gene_ids")) + _split_ids(row.get("lncRNA_gene_ids"))
        own[locus_id] = {i.split(".")[0] for i in ids}
    for name in control_genes:
        own["ctrl:" + name] = set(ctrl.loc[ctrl["gene_name"] == name, "gene_id"].str.split(".").str[0])
    g = gencode_exons[["chrom", "start0", "end1"]].assign(gene_id=gencode_exons["gene_id"].str.split(".").str[0])
    g = g.drop_duplicates()
    g_by_chrom = dict(tuple(g.groupby("chrom")))
    regions = []
    for (fid, chrom), sub in exons.groupby(["feature_id", "chrom"]):
        blockers = g_by_chrom.get(chrom, g.iloc[:0])
        blockers = blockers[(blockers["start0"] < sub["end1"].max()) & (blockers["end1"] > sub["start0"].min())]
        blockers = blockers[~blockers["gene_id"].isin(own.get(fid, set()))]
        regions.append(subtract_intervals(sub, blockers))
    regions = merge_intervals(pd.concat(regions, ignore_index=True))

    exonic = exons.assign(n=exons["end1"] - exons["start0"]).groupby("feature_id")["n"].sum()
    quant = regions.assign(n=regions["end1"] - regions["start0"]).groupby("feature_id")["n"].sum()
    span = exons.groupby("feature_id").agg(chrom=("chrom", "first"), start0=("start0", "min"), end1=("end1", "max"))
    features = span.join(pd.DataFrame({"exonic_bases": exonic})).join(quant.rename("quant_bases")).fillna({"quant_bases": 0})
    features["quant_bases"] = features["quant_bases"].astype(int)
    features["quant_fraction"] = features["quant_bases"] / features["exonic_bases"]
    features["quantifiable"] = features["quant_bases"] >= min_quant_bases
    return regions, features


# ---------------------------------------------------------------- samples

def _metadata_url(source, project, table):
    sub = project[-2:]
    return f"{RECOUNT3}{source}/metadata/{sub}/{project}/{source}.{table}.{project}.MD.gz"


def _read_url_tsv(url):
    with urllib.request.urlopen(url, timeout=300) as r:
        return pd.read_csv(gzip.GzipFile(fileobj=r), sep="\t", low_memory=False)


def bigwig_url(source, project, external_id):
    # the sub-directory is named after the last two characters of the id, upper-case and without a version
    # suffix (GTEx '.1')
    stem = re.sub(r"\.\d+$", "", external_id)
    return (f"{RECOUNT3}{source}/base_sums/{project[-2:]}/{project}/{stem[-2:].upper()}/"
            f"{source}.base_sums.{project}_{external_id}.ALL.bw")


def gene_sums_url(source, project, annotation="G026"):
    return f"{RECOUNT3}{source}/gene_sums/{project[-2:]}/{project}/{source}.gene_sums.{project}.{annotation}.gz"


def read_gene_sums(path, source, samples):
    """recount3 gene sums (GENCODE v26) of the given samples, divided by their AUC; index = unversioned gene id,
    columns = sample_id."""
    g = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    g.columns = source + ":" + g.columns.astype(str)
    auc = samples.set_index("sample_id")["auc"]
    g = g[[c for c in g.columns if c in auc.index]]
    g = g / auc[g.columns] * 1e9
    g.index = g.index.str.split(".").str[0]
    return g[~g.index.duplicated()]


def manifest(source, project, keep_columns=()):
    """
    Samples of a recount3 project (source tcga / gtex / sra) with bigWig URL, total coverage (AUC) and the
    requested phenotype columns of the source-specific metadata table.
    """
    qc = _read_url_tsv(_metadata_url(source, project, "recount_qc"))
    pheno = _read_url_tsv(_metadata_url(source, project, source))
    m = qc[["rail_id", "external_id", "study", "bc_auc.all_reads_all_bases"]].rename(
        columns={"bc_auc.all_reads_all_bases": "auc"})
    keep = [c for c in keep_columns if c in pheno.columns]
    on = "rail_id" if "rail_id" in pheno.columns else "external_id"
    m = m.merge(pheno[[on] + [c for c in keep if c != on]], on=on, how="left")
    m.insert(0, "source", source)
    m.insert(1, "project", project)
    m["sample_id"] = source + ":" + m["external_id"].astype(str)
    m["url"] = [bigwig_url(source, project, str(e)) for e in m["external_id"]]
    return m


# ---------------------------------------------------------------- quantification

def sum_regions(bw, regions):
    """Summed coverage per feature over its regions in an open pyBigWig file (missing chromosomes -> 0)."""
    chroms = bw.chroms()
    sums = np.zeros(len(regions))
    for i, (c, s, e) in enumerate(zip(regions["chrom"], regions["start0"], regions["end1"])):
        if c in chroms:
            v = bw.stats(c, int(s), int(min(e, chroms[c])), type="sum", exact=True)[0]
            sums[i] = v or 0.0
    return pd.Series(sums, index=regions.index).groupby(regions["feature_id"]).sum()


def _download(url, dest):
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 22)
    os.replace(tmp, dest)


def run(samples, regions, out_dir, worker="0", scratch=None):
    """
    Stream the bigWigs of `samples` (manifest rows); append one row of sums per sample, skip done ones. The regions
    are kept in out_dir/regions.tsv; a later run on changed regions raises (expression.read_sums checks the same).
    Samples that fail (download, bigWig) are reported and the run exits non-zero after the others, so that a rerun
    retries them.
    """
    import pyBigWig
    from .expression import REGION_COLUMNS, changed_features
    os.makedirs(out_dir, exist_ok=True)
    saved = os.path.join(out_dir, "regions.tsv")
    if os.path.exists(saved):
        changed = changed_features(pd.read_csv(saved, sep="\t"), regions)
        if changed:
            raise ValueError(f"{out_dir}: regions of {len(changed)} features differ from those of the existing sums")
    else:
        regions[REGION_COLUMNS].to_csv(saved, sep="\t", index=False)
    features = sorted(regions["feature_id"].unique())
    path = os.path.join(out_dir, f"sums.{worker}.tsv")
    done = set()
    if os.path.exists(path):
        done = set(pd.read_csv(path, sep="\t", usecols=["sample_id"])["sample_id"])
    else:
        with open(path, "w") as fh:
            fh.write("\t".join(["sample_id"] + features) + "\n")
    tmpdir = tempfile.mkdtemp(dir=scratch)
    failed = []
    try:
        for row in samples.itertuples(index=False):
            if row.sample_id in done:
                continue
            local = os.path.join(tmpdir, "x.bw")
            try:
                _download(row.url, local)
                bw = pyBigWig.open(local)
                sums = sum_regions(bw, regions).reindex(features).fillna(0.0)
                bw.close()
            except Exception as exc:          # one bad sample must not stop the stream
                print(f"failed {row.sample_id}: {exc}", flush=True)
                failed.append(row.sample_id)
                continue
            finally:
                if os.path.exists(local):
                    os.remove(local)
            with open(path, "a") as fh:
                fh.write("\t".join([row.sample_id] + [f"{v:.6g}" for v in sums]) + "\n")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    if failed:
        raise SystemExit(f"{len(failed)} samples failed ({', '.join(failed[:5])}); rerun to retry them")
    return path


def normalize(sums, features, samples):
    """ncov = summed coverage / quantifiable positions / sample AUC * 1e9 (features x samples)."""
    s = sums.set_index("sample_id").T
    n = features["quant_bases"].reindex(s.index).replace(0, np.nan)
    auc = samples.set_index("sample_id")["auc"].reindex(s.columns)
    return s.div(n, axis=0).div(auc, axis=1) * 1e9
