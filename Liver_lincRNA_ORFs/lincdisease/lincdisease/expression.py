"""
Expression of the quantified features (LOGIC.md, Step 2): the coverage sums written by quant.run, checked against
the regions they were computed on, turned into ncov (quant.normalize) and scaled by size factors.

Sums and regions. Every sums directory (<sample set>/sums, <sample set>/sums_<name> for features added later)
keeps a copy of the regions it was computed on (regions.tsv, written by quant.run). read_sums raises if the regions
of a feature have changed since, or if a feature of the current regions has no sums: ncov would otherwise be
computed from coverage of other positions, or be missing without notice.

Sample scaling. ncov divides by the total base coverage (AUC) of a sample. In liver a few hepatocyte transcripts
(ALB, APOA1, ...) take a large and variable share of it, which falls in tumours and injured liver and so lifts
every other transcript. ncov is therefore divided by median-of-ratios size factors of the protein-coding recount3
gene sums of the same samples (size_factors), computed over all samples that are compared with each other.
"""
import glob
import os

import numpy as np
import pandas as pd

from . import quant

REGION_COLUMNS = ["feature_id", "chrom", "start0", "end1"]
# size factors: genes detected in at least this fraction of the samples, the most abundant 1% left out (they are
# the composition-driven transcripts themselves)
SF_MIN_DETECTED = 0.9
SF_TOP_EXCLUDED = 0.01


def size_factors(genes):
    """
    Median-of-ratios size factors (geometric mean 1) of gene sums scaled by total coverage (genes x samples,
    protein-coding). Equivalent to DESeq2's estimateSizeFactors on the kept genes: per sample, the median over genes
    of log(y_gj) - mean_j log(y_gj), exponentiated; zeros are left out per gene (log of 0 is missing).
    """
    g = genes[(genes > 0).mean(axis=1) >= SF_MIN_DETECTED]
    g = g.drop(g.mean(axis=1).nlargest(max(1, int(len(g) * SF_TOP_EXCLUDED))).index)
    lg = np.log(g.where(g > 0))
    sf = np.exp(lg.sub(lg.mean(axis=1), axis=0).median())
    return sf / np.exp(np.log(sf).mean())


def _region_set(regions):
    """feature_id -> sorted tuple of (chrom, start0, end1) of its regions."""
    r = regions[REGION_COLUMNS].sort_values(REGION_COLUMNS)
    return {f: tuple(map(tuple, g[["chrom", "start0", "end1"]].to_numpy())) for f, g in r.groupby("feature_id")}


def changed_features(saved, current):
    """Features present in both region tables whose regions differ."""
    a, b = _region_set(saved), _region_set(current)
    return sorted(f for f in a.keys() & b.keys() if a[f] != b[f])


def read_sums(sample_dir, regions):
    """
    Coverage sums of <sample_dir>/sums*/ joined by sample, for the features of `regions` (current quant/regions.tsv).
    Raises if a sums directory has no regions copy, if regions changed since a directory was computed, or if a
    feature has no sums in some or all samples.
    """
    parts = []
    for sd in sorted(glob.glob(os.path.join(sample_dir, "sums*"))):
        files = sorted(glob.glob(os.path.join(sd, "sums.*.tsv")))
        if not files:
            continue
        saved = os.path.join(sd, "regions.tsv")
        if not os.path.exists(saved):
            raise ValueError(f"{sd} has no regions.tsv (the regions its sums were computed on)")
        changed = changed_features(pd.read_csv(saved, sep="\t"), regions)
        if changed:
            raise ValueError(f"{sd}: the regions of {len(changed)} features changed since the sums were computed "
                             f"(e.g. {changed[:3]}); delete the directory to quantify again")
        parts.append(pd.concat([pd.read_csv(f, sep="\t") for f in files]).drop_duplicates("sample_id")
                     .set_index("sample_id"))
    if not parts:
        raise FileNotFoundError(f"no coverage sums in {sample_dir}/sums*/")
    sums = parts[0]
    for extra in parts[1:]:
        sums = sums.join(extra[extra.columns.difference(sums.columns)], how="left")
    wanted = sorted(regions["feature_id"].unique())
    missing = [f for f in wanted if f not in sums.columns]
    if missing:
        raise ValueError(f"{sample_dir}: no sums for {len(missing)} features ({missing[:3]}); quantify them with "
                         f"`quant --sample-set {os.path.basename(sample_dir)} --features ... --sums-dir sums_<name>`")
    sums = sums.copy()[wanted]                      # copy: the joins leave a fragmented frame
    gaps = sums.isna()
    if gaps.any().any():     # a sample missing from a sums_<name> directory, or a feature no longer quantified
        samples, feats = gaps.any(axis=1), gaps.any()
        raise ValueError(f"{sample_dir}: no sums for {int(samples.sum())} samples x {int(feats.sum())} features "
                         f"(e.g. {list(samples[samples].index[:3])} x {list(feats[feats].index[:3])}); quantify them "
                         f"with `quant --sample-set {os.path.basename(sample_dir)} --features ...`")
    return sums.reset_index()


def quantifiable(feats):
    """Features measurable in unstranded public RNA-seq (Step 2); the others are not tested."""
    return feats.index[feats["quantifiable"].fillna(False).astype(bool)]


def load(out_dir, sample_set="quant"):
    """AUC-scaled ncov (features x samples) of a sample set with its sample table and quant/features.tsv."""
    regions = pd.read_csv(os.path.join(out_dir, "quant", "regions.tsv"), sep="\t")
    feats = pd.read_csv(os.path.join(out_dir, "quant", "features.tsv"), sep="\t", index_col=0)
    d = os.path.join(out_dir, sample_set)
    sums = read_sums(d, regions)
    samples = pd.read_csv(os.path.join(d, "samples.tsv"), sep="\t", low_memory=False)
    samples = samples[samples["sample_id"].isin(sums["sample_id"])]
    return quant.normalize(sums, feats, samples), samples, feats


def protein_coding_ids(gtf_path):
    """Unversioned gene ids of the protein-coding genes of a GENCODE GTF."""
    g = quant.read_gtf_exons(gtf_path).drop_duplicates("gene_id")
    return set(g.loc[g["gene_type"] == "protein_coding", "gene_id"].str.split(".").str[0])


def scaled(ncov, sf):
    """ncov divided by the size factor of each sample (columns without one are dropped)."""
    cols = ncov.columns.intersection(sf.index)
    return ncov[cols].div(sf[cols], axis=1)


def read_size_factors(out_dir, sample_set="quant"):
    return pd.read_csv(os.path.join(out_dir, sample_set, "size_factors.tsv"), sep="\t", index_col=0)["size_factor"]
