"""
Mechanism classes (LOGIC.md, Step 4): per-ORF features (with locus-level features copied to the locus' ORFs),
percentiles against all ORFs of the run, class scores, relative shares and dominant class.

Features computed here (higher = more support for the class; '-' marks features ranked on the negated value):
  M1  translation (pooled log10 Bayes factor), fraction of the transcript's P-sites in the ORF, specific
      peptides (proteomics peptides absent from canonical proteins + lincORF-specific HLA-I peptides),
      PhyloCSF weighted by copy-attribution confidence, ORF length
  M2  NMD substrate (stop codon > 50 nt upstream of the last exon junction), number of other translated ORFs
      (evidence moderate/strong) on the locus
  M3  divergent / antisense position next to a protein-coding gene, low polysome association (-MRL percentile),
      co-expression with the nearest protein-coding gene (TSS distance <= 100 kb) across TCGA-LIHC and GTEx
      liver: Spearman rho minus the median rho of protein-coding gene pairs in the same TSS-distance bin
  M4  copy number (log2 TPM HepG2), no cis position (1 - M3 position)
  M5  ribosome association without one dominant ORF (MRL percentile x (1 - largest ORF share of the locus))
  M6  snoRNA / scaRNA / miRNA genes (GENCODE + UCSC wgRna) inside the locus span on the same strand
ENCODE HepG2 features (encode.py): nuclear enrichment (M3, M4) and its negative (M5). UPF1 knock-down (M2),
eCLIP, chromatin at the TSS, ORF vs non-ORF exon conservation, sequence features (TM, signal peptide, Pfam,
disorder) and CRISPR screens enter as missing until their inputs are added; scoring.composite shrinks them to 0.5.
"""
import numpy as np
import pandas as pd

from . import quant, scoring

CLASSES = {
    "M1": ["m1_translation", "m1_orf_share", "m1_peptides", "m1_phylocsf_weighted", "m1_length",
           "m1_orf_vs_exon_conservation", "m1_sequence_features", "m1_crispr_screen"],
    "M2": ["m2_nmd", "m2_other_translated_orfs", "m2_upf1_kd"],
    "M3": ["m3_neighbour_position", "m3_low_polysome", "m3_coexpression", "m3_nuclear", "m3_enhancer_chromatin"],
    "M4": ["m4_copy_number", "m4_no_cis", "m4_nuclear", "m4_eclip_nuclear"],
    "M5": ["m5_ribosome_no_dominant_orf", "m5_cytoplasmic", "m5_eclip_cytoplasmic", "m5_mirna_sites"],
    "M6": ["m6_small_rna_host"],
}
DOMINANCE_MARGIN = 0.05
DIVERGENT_DISTANCE = 5000
NEIGHBOUR_MAX_DISTANCE = 100_000
DISTANCE_BINS = (0, 5_000, 20_000, 50_000, 100_000)
SMALL_RNA_TYPES = {"snoRNA", "scaRNA", "miRNA"}


# ---------------------------------------------------------------- peptides

def specific_peptides(proteomics_peptides, proteomics_entries, immunopeptides, immuno_entries):
    """
    Per orf_id: n_specific_peptides (proteomics peptides not in canonical proteins/contaminants), n_hla_specific
    (lincORF-specific HLA-I peptides) and copy_specific (a peptide of either kind that maps to this ORF only).
    Each search has its own entry table; an entry may stand for several ORFs with identical sequences.
    """
    def explode(df, peptide_col, entries):
        d = df[[peptide_col, "lincORF_entries"]].dropna()
        d = d.assign(entry=d["lincORF_entries"].str.split(r"[;,]")).explode("entry")
        d["entry"] = d["entry"].str.strip()
        d = d.merge(entries[["entry", "orf_id"]], on="entry")
        d["n_orfs"] = d.groupby(peptide_col)["orf_id"].transform("nunique")
        return d

    p = explode(proteomics_peptides[~proteomics_peptides["in_canonical_or_contaminant"].astype(str).eq("True")],
                "peptide", proteomics_entries)
    h = explode(immunopeptides[immunopeptides["lincORF_specific"].astype(str).eq("True")], "Peptide",
                immuno_entries)
    out = pd.DataFrame({
        "n_specific_peptides": p.groupby("orf_id")["peptide"].nunique(),
        "n_hla_specific": h.groupby("orf_id")["Peptide"].nunique(),
    }).fillna(0)
    single = set(p.loc[p["n_orfs"] == 1, "orf_id"]) | set(h.loc[h["n_orfs"] == 1, "orf_id"])
    out["copy_specific"] = out.index.isin(single)
    return out


# ---------------------------------------------------------------- structure

def stop_to_last_junction(orfs, tx_exons):
    """Exonic distance (nt) from the ORF stop codon end to the last exon junction of its transcript (0 if none)."""
    ex = {t: g.sort_values("start0") for t, g in tx_exons.groupby("transcript_id")}
    out = pd.Series(np.nan, index=orfs.index)
    for idx, r in orfs[["transcript_id", "orf_start0", "orf_end1", "strand"]].iterrows():
        g = ex.get(r["transcript_id"])
        if g is None:
            continue
        if len(g) == 1:
            out[idx] = 0
            continue
        s, e = g["start0"].to_numpy(), g["end1"].to_numpy()
        if r["strand"] == "+":
            stop, junction = r["orf_end1"], s[-1]                   # start of the last exon
            d = np.clip(np.minimum(e, junction) - np.maximum(s, stop), 0, None).sum()
        else:
            stop, junction = r["orf_start0"], e[0]                  # end of the last exon (leftmost)
            d = np.clip(np.minimum(e, stop) - np.maximum(s, junction), 0, None).sum()
        out[idx] = d
    return out


def gene_spans(df, id_col):
    """Span, strand and TSS per `id_col` of an exon/transcript table."""
    g = df.groupby(id_col).agg(chrom=("chrom", "first"), start0=("start0", "min"), end1=("end1", "max"),
                               strand=("strand", "first"))
    g["tss"] = np.where(g["strand"] == "+", g["start0"], g["end1"])
    return g


def neighbour_position(loci_tx, gencode_exons, distance=DIVERGENT_DISTANCE):
    """
    Per locus: 1 if it overlaps exons of a protein-coding gene on the opposite strand (antisense) or its TSS lies
    within `distance` of the TSS of a protein-coding gene on the opposite strand, upstream (divergent).
    `loci_tx`: lincRNA isoform table (gene_id, chrom, start0, end1, strand).
    """
    pc = gencode_exons[gencode_exons["gene_type"] == "protein_coding"]
    genes = gene_spans(pc, "gene_id")
    span = gene_spans(loci_tx, "gene_id")
    pc_by = {(c, s): g for (c, s), g in pc.groupby(["chrom", "strand"])}
    genes_by = {(c, s): g for (c, s), g in genes.groupby(["chrom", "strand"])}
    out = pd.Series(0.0, index=span.index)
    for lid, r in span.iterrows():
        opp = "-" if r["strand"] == "+" else "+"
        ex = pc_by.get((r["chrom"], opp))
        if ex is not None and ((ex["start0"] < r["end1"]) & (ex["end1"] > r["start0"])).any():
            out[lid] = 1.0
            continue
        g = genes_by.get((r["chrom"], opp))
        if g is None:
            continue
        # divergent: the coding gene's TSS lies upstream of the lincRNA TSS (head-to-head)
        d = (r["tss"] - g["tss"]) if r["strand"] == "+" else (g["tss"] - r["tss"])
        if ((d >= 0) & (d <= distance)).any():
            out[lid] = 1.0
    return out


def nearest_coding(query, genes, max_distance=NEIGHBOUR_MAX_DISTANCE, exclude_self=False):
    """
    For each row of `query` (index id; chrom, tss) the protein-coding gene of `genes` (same layout) with the
    nearest TSS, on either strand, within `max_distance`: columns neighbour (unversioned gene id), distance.
    """
    out = pd.DataFrame({"neighbour": pd.Series(np.nan, index=query.index, dtype=object),
                        "distance": np.nan})
    by = {c: g.sort_values("tss") for c, g in genes.groupby("chrom")}
    for chrom, q in query.groupby("chrom"):
        g = by.get(chrom)
        if g is None:
            continue
        pos, ids = g["tss"].to_numpy(), g.index.to_numpy()
        for qid, t in q["tss"].items():
            k = np.searchsorted(pos, t)
            cand = [j for j in (k - 2, k - 1, k, k + 1) if 0 <= j < len(pos) and not (exclude_self and ids[j] == qid)]
            j = min(cand, key=lambda j: abs(pos[j] - t), default=None)
            if j is not None and abs(pos[j] - t) <= max_distance:
                out.loc[qid, ["neighbour", "distance"]] = [ids[j], abs(pos[j] - t)]
    return out


def coexpression(locus_expr, gene_expr, neighbours, null_pairs, bins=DISTANCE_BINS):
    """
    Spearman rho between each locus and its neighbour gene over the shared samples, minus the median rho of
    `null_pairs` (gene -> neighbour gene, distance) in the same distance bin. Returns rho, null and excess.
    """
    ranks_l = locus_expr.rank(axis=1)
    ranks_g = gene_expr.rank(axis=1)

    def rho(a, b):
        a, b = a - a.mean(), b - b.mean()
        den = np.sqrt((a ** 2).sum() * (b ** 2).sum())
        return float((a * b).sum() / den) if den > 0 else np.nan

    def bin_of(d):
        return np.digitize(d, bins[1:-1])

    nul = null_pairs.dropna(subset=["neighbour"])
    nul = nul[nul.index.isin(ranks_g.index) & nul["neighbour"].isin(ranks_g.index)]
    nul_rho = [rho(ranks_g.loc[a].to_numpy(), ranks_g.loc[b].to_numpy()) for a, b in nul["neighbour"].items()]
    null_median = pd.Series(nul_rho, index=nul.index).groupby(bin_of(nul["distance"].to_numpy())).median()

    out = pd.DataFrame(index=locus_expr.index, columns=["neighbour", "distance", "coexpr_rho", "coexpr_null"])
    for lid, r in neighbours.dropna(subset=["neighbour"]).iterrows():
        if lid not in ranks_l.index or r["neighbour"] not in ranks_g.index:
            continue
        out.loc[lid] = [r["neighbour"], r["distance"],
                        rho(ranks_l.loc[lid].to_numpy(), ranks_g.loc[r["neighbour"]].to_numpy()),
                        null_median.get(bin_of(r["distance"]), np.nan)]
    out = out.astype({"distance": float, "coexpr_rho": float, "coexpr_null": float})
    out["coexpr_excess"] = out["coexpr_rho"] - out["coexpr_null"]
    return out


def small_rna_genes(gencode_exons, wgrna):
    """
    snoRNA / scaRNA / miRNA genes of GENCODE and of UCSC wgRna (snoRNABase + miRBase; GENCODE v42 lacks e.g. the
    GAS5 snoRNAs), merged so that a gene present in both counts once: chrom/start0/end1/strand.
    """
    g = gencode_exons[gencode_exons["gene_type"].isin(SMALL_RNA_TYPES)]
    g = g.groupby("gene_id").agg(chrom=("chrom", "first"), start0=("start0", "min"), end1=("end1", "max"),
                                 strand=("strand", "first"))
    both = pd.concat([g[["chrom", "start0", "end1", "strand"]], wgrna[["chrom", "start0", "end1", "strand"]]])
    merged = quant.merge_intervals(both.rename(columns={"strand": "feature_id"}))
    return merged.rename(columns={"feature_id": "strand"})


def small_rna_hosts(loci_tx, small_rnas):
    """Per locus: number of small-RNA genes (small_rna_genes) within its span on the same strand."""
    span = gene_spans(loci_tx, "gene_id")
    by = {(c, s): g for (c, s), g in small_rnas.groupby(["chrom", "strand"])}
    out = pd.Series(0.0, index=span.index)
    for lid, r in span.iterrows():
        g = by.get((r["chrom"], r["strand"]))
        if g is not None:
            out[lid] = float(((g["start0"] >= r["start0"]) & (g["end1"] <= r["end1"])).sum())
    return out


# ---------------------------------------------------------------- features and scores

def features(orfs, qc_table, peptides, stop_dist, position, small_rna, coexpr=None, encode=None):
    """
    Raw feature values per ORF (index of `orfs`). `coexpr`: locus -> co-expression excess (coexpression);
    `encode`: locus table of encode.features.
    """
    f = pd.DataFrame(index=orfs.index)
    num = lambda c: pd.to_numeric(orfs[c], errors="coerce")
    pep = peptides.reindex(orfs["orf_id"]).fillna(0)
    f["m1_translation"] = num("bayes_log10BF_pooled")
    f["m1_orf_share"] = num("frac_tx_psites_in_orf")
    f["m1_peptides"] = (pep["n_specific_peptides"] + pep["n_hla_specific"]).to_numpy()
    f["m1_phylocsf"] = num("PhyloCSF_Median")
    f["copy_attribution"] = qc_table["copy_attribution"].to_numpy()
    f["m1_length"] = num("n_aa")
    f["m2_nmd"] = (stop_dist > 50).astype(float).where(stop_dist.notna())
    translated = orfs["ribo_evidence"].isin(["moderate", "strong"])
    per_locus = translated.groupby(orfs["locus_id"]).transform("sum")
    f["m2_other_translated_orfs"] = (per_locus - translated.astype(int)).astype(float)
    f["m3_neighbour_position"] = orfs["locus_id"].map(position).to_numpy()
    f["m3_low_polysome"] = -num("MRL_pct_vs_mRNA")
    f["m3_coexpression"] = orfs["locus_id"].map(coexpr).to_numpy() if coexpr is not None else np.nan
    enc = encode if encode is not None else pd.DataFrame(columns=["nuclear_log2ratio"])
    nuclear = orfs["locus_id"].map(enc["nuclear_log2ratio"]).astype(float).to_numpy()
    # m2_upf1_kd stays missing: the ENCODE HepG2 UPF1 shRNA knock-down does not inhibit NMD (encode.py)
    f["m3_nuclear"] = nuclear
    f["m4_copy_number"] = num("log2TPM_HepG2")
    f["m4_no_cis"] = 1.0 - f["m3_neighbour_position"]
    f["m4_nuclear"] = nuclear
    max_share = f["m1_orf_share"].groupby(orfs["locus_id"]).transform("max")
    f["m5_ribosome_no_dominant_orf"] = num("MRL_pct_vs_mRNA") * (1.0 - max_share)
    f["m5_cytoplasmic"] = -nuclear
    f["m6_small_rna_host"] = orfs["locus_id"].map(small_rna).to_numpy()
    return f


def class_table(feats, is_background, margin=DOMINANCE_MARGIN):
    """Percentiles, class scores, coverage, relative shares and dominant class."""
    cols = [c for c in feats.columns if c not in ("copy_attribution", "m1_phylocsf")]
    pct = scoring.component_percentiles(feats, feats[is_background], cols)
    # PhyloCSF counts in proportion to the confidence that the alignment belongs to this copy
    p_cons = scoring.percentile_vs_background(feats["m1_phylocsf"], feats.loc[is_background, "m1_phylocsf"])
    pct["m1_phylocsf_weighted"] = 0.5 + (p_cons - 0.5) * feats["copy_attribution"].to_numpy()
    scores, coverage = scoring.class_scores(pct, CLASSES)
    shares = scoring.relative_shares(scores)
    out = pd.concat([pct.add_prefix("pct_"), scores.add_suffix("_score"), coverage.add_suffix("_coverage"),
                     shares.add_suffix("_share")], axis=1)
    out["dominant_class"] = scoring.dominant_class(shares, margin)
    return out
