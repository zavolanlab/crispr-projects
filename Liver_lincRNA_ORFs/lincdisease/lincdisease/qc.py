"""
Quality-check layer (LOGIC.md, Step 1): flags, copy-attribution confidence and confidence tier per ORF.
Flags never remove ORFs and never change scores; copy-attribution confidence weights the M1 features a paralogous
copy can supply (translation, P-site share, sequence features, conservation: mechanism.COPY_WEIGHTED).

Copy-attribution confidence (0..1) applies to ORFs in segmental duplications / falsely duplicated GRCh38
regions or with homology to annotated proteins (other ORFs get 1). It is the mean of the available parts:
  expression   fraction of P-sites from uniquely mapped reads
  peptides     1 if a peptide covers a residue that differs from the paralogs (missing otherwise: absence of
               peptides is not evidence against the copy)
  age          (identity_max - identity) / (identity_max - identity_min), clipped to 0..1: very recent copies
               (identity >= identity_max) leave the multi-species alignment with the parent's orthologs
"""
import numpy as np
import pandas as pd

TIERS = ("clean", "caution", "expressed_paralog", "likely_artefact")


def _nonempty(series):
    return series.fillna("").astype(str).str.strip().ne("")


def _as_bool(series):
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def orf_exon_blocks(orfs):
    """One row per coding exon block of each ORF (index = ORF index; chrom, orf_start0, orf_end1), from the lincorfz
    `orf_blocks` ('start0-end1;...'); ORFs without blocks keep their span."""
    if "orf_blocks" not in orfs:
        return orfs[["chrom", "orf_start0", "orf_end1"]]
    rows = []
    for idx, chrom, blocks, s0, e1 in zip(orfs.index, orfs["chrom"], orfs["orf_blocks"], orfs["orf_start0"],
                                          orfs["orf_end1"]):
        parts = [b.split("-") for b in str(blocks).split(";") if "-" in b] if isinstance(blocks, str) else []
        rows += [(idx, chrom, int(a), int(b)) for a, b in parts] or [(idx, chrom, s0, e1)]
    return pd.DataFrame(rows, columns=["idx", "chrom", "orf_start0", "orf_end1"]).set_index("idx")


def interval_overlap_max(orfs, regions, value_col):
    """
    For each ORF the maximum `value_col` of regions (chrom, start0, end1) overlapping its coding exons (introns of
    multi-exon ORFs do not count), NaN without overlap. Strand is ignored (duplications are copy-level).
    """
    out = pd.Series(np.nan, index=orfs.index, dtype=float)
    blocks = orf_exon_blocks(orfs)
    for chrom, sub in blocks.groupby("chrom"):
        reg = regions[regions["chrom"] == chrom]
        if reg.empty:
            continue
        rs = reg["start0"].to_numpy()
        re_ = reg["end1"].to_numpy()
        rv = reg[value_col].to_numpy(dtype=float)
        order = np.argsort(rs)
        rs, re_, rv = rs[order], re_[order], rv[order]
        # regions are long and may nest, so scan all regions starting before the ORF end
        for idx, s, e in zip(sub.index, sub["orf_start0"].to_numpy(), sub["orf_end1"].to_numpy()):
            k = np.searchsorted(rs, e, side="left")
            hit = re_[:k] > s
            if hit.any():
                out[idx] = np.nanmax([out[idx], rv[:k][hit].max()])
    return out


def qc_flags(orfs, segdups=None, false_dups=None, homology=None, quantifiable=None):
    """
    Collect the QC flags of Step 1.
      segdups       DataFrame chrom/start0/end1/fracMatch (UCSC genomicSuperDups)
      false_dups    DataFrame chrom/start0/end1 (GRCh38 regions falsely duplicated relative to T2T-CHM13)
      homology      DataFrame orf_id/parent/pident/qcov (best DIAMOND hit vs annotated proteins)
      quantifiable  Series orf_id -> bool (enough exonic positions not overlapping the sense gene; Step 2)
    """
    f = pd.DataFrame(index=orfs.index)
    f["orf_id"] = orfs["orf_id"]
    has = lambda table: table is not None and len(table) > 0
    f["segdup_fracMatch"] = interval_overlap_max(orfs, segdups, "fracMatch") if has(segdups) else np.nan
    f["grch38_false_dup"] = (interval_overlap_max(orfs, false_dups.assign(one=1.0), "one").eq(1.0)
                             if has(false_dups) else False)
    f["in_duplication"] = f["segdup_fracMatch"].notna() | f["grch38_false_dup"]

    if has(homology):
        h = homology.set_index("orf_id")
        f["homolog_parent"] = orfs["orf_id"].map(h["parent"])
        f["homolog_pident"] = orfs["orf_id"].map(h["pident"]).astype(float)
        f["homolog_qcov"] = orfs["orf_id"].map(h["qcov"]).astype(float)
    else:
        f[["homolog_parent", "homolog_pident", "homolog_qcov"]] = np.nan
    f["protein_homolog"] = f["homolog_parent"].notna()

    f["antisense_to_coding"] = (_nonempty(orfs["orf_antisense_to_CDS"])
                                | _nonempty(orfs["antisense_coding_exon_overlap"]))
    f["small_RNA_overlap"] = _nonempty(orfs["orf_overlaps_small_RNA"])
    f["multimapper_dominated"] = _as_bool(orfs["resolved_multimapper_dominated"])
    f["coding_gene_extension"] = (_nonempty(orfs["possible_3UTR_extension_of"])
                                  | _nonempty(orfs["possible_5prime_extension_of"]))
    f["tcga_quantifiable"] = (orfs["orf_id"].map(quantifiable).fillna(False).astype(bool) if quantifiable is not None
                              else np.nan)
    return f


def copy_attribution(orfs, flags, copy_peptides=None, identity_max=0.98, identity_min=0.90):
    """
    Copy-attribution confidence per ORF (see module docstring). Also returns the parts.
      copy_peptides  Series orf_id -> bool, peptide on a paralog-discriminating residue
    """
    parts = pd.DataFrame(index=orfs.index)
    parts["attr_expression"] = pd.to_numeric(orfs["frac_psites_unique"], errors="coerce").clip(0, 1)
    if copy_peptides is not None:
        pep = orfs["orf_id"].map(copy_peptides).fillna(False).astype(bool)
        parts["attr_peptides"] = np.where(pep, 1.0, np.nan)
    else:
        parts["attr_peptides"] = np.nan
    identity = flags["segdup_fracMatch"].fillna(flags["homolog_pident"] / 100.0)
    parts["attr_age"] = ((identity_max - identity) / (identity_max - identity_min)).clip(0, 1)

    affected = flags["in_duplication"] | flags["protein_homolog"]
    conf = parts.mean(axis=1, skipna=True)
    parts["copy_attribution"] = np.where(affected, conf.fillna(0.5), 1.0)
    return parts


def confidence_tier(flags, parts, min_unique=0.5):
    """
    likely_artefact    duplicated/homologous, and the translation evidence cannot be attributed to the copy:
                       < `min_unique` unique P-sites (or multimapper-dominated) and no copy-specific peptide
    expressed_paralog  duplicated/homologous otherwise
    caution            antisense to coding exons, small-RNA overlap, multimapper-dominated, coding-gene
                       extension or not quantifiable in TCGA/GTEx
    clean              otherwise
    """
    affected = flags["in_duplication"] | flags["protein_homolog"]
    no_unique = parts["attr_expression"].fillna(0).lt(min_unique) | flags["multimapper_dominated"]
    no_pep = parts["attr_peptides"].isna()
    artefact = affected & no_unique & no_pep
    not_quant = flags["tcga_quantifiable"].eq(False)
    caution = (flags["antisense_to_coding"] | flags["small_RNA_overlap"] | flags["multimapper_dominated"]
               | flags["coding_gene_extension"] | not_quant)
    tier = pd.Series("clean", index=flags.index, dtype=object)
    tier[caution] = "caution"
    tier[affected] = "expressed_paralog"
    tier[artefact] = "likely_artefact"
    return tier


def run_qc(orfs, **kwargs):
    """Flags, copy-attribution parts and confidence tier in one table (index of `orfs`)."""
    flag_kw = {k: kwargs[k] for k in ("segdups", "false_dups", "homology", "quantifiable") if k in kwargs}
    attr_kw = {k: kwargs[k] for k in ("copy_peptides",) if k in kwargs}
    flags = qc_flags(orfs, **flag_kw)
    parts = copy_attribution(orfs, flags, **attr_kw)
    out = pd.concat([flags, parts], axis=1)
    out["confidence_tier"] = confidence_tier(flags, parts)
    return out
