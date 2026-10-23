#!/usr/bin/env python
"""
Coverage-weighted mean ribosome load (MRL) of lincRNA isoforms, ORFs and reference mRNAs.

For a feature (set of exonic intervals on one strand), C_k is the summed replicate-mean
polysome coverage of the fraction with k ribosomes over the feature, and

    MRL = sum_k k * C_k / sum_k C_k       (reported if sum_k C_k / length >= --min-mean-coverage)

Polysome profiling sediments whole RNA molecules, so the isoform-level MRL (all exons) is the
primary estimate of how many ribosomes a lincRNA molecule carries; the ORF-level MRL uses only
the ORF blocks and is reported for comparison. Ensembl_canonical protein-coding transcripts
serve as the reference distribution (MRL_pct_vs_mRNA = percentile among expressed mRNAs).

Outputs: polysome.transcripts.tsv, polysome.orfs.tsv, polysome.reference_mRNAs.tsv, polysome_mrl.pdf
"""
import argparse
import multiprocessing as mp

import matplotlib
import numpy as np
import pandas as pd
import pyBigWig

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from lincorfz_utils import bigwig_chrom_name, read_chromosomes, read_gtf  # noqa: E402

G = {}


def init_worker(track_table):
    G["bw"] = {(int(r.n_ribosomes), r.strand): pyBigWig.open(r.path) for r in track_table.itertuples()}
    G["ks"] = sorted({k for k, _ in G["bw"]})


def feature_sums(rec):
    fid, chrom, strand, blocks = rec
    sums = []
    for k in G["ks"]:
        bw = G["bw"][(k, strand)]
        c = bigwig_chrom_name(bw, chrom)
        tot = 0.0
        if c is not None:
            for s, e in blocks:
                v = bw.stats(c, int(s), int(e), type="sum", exact=True)[0]
                tot += v if v is not None and np.isfinite(v) else 0.0
        sums.append(tot)
    length = int(sum(e - s for s, e in blocks))
    return fid, length, sums


def summarize(results, ks, min_mean_cov, prefix=""):
    rows = []
    for fid, length, sums in results:
        sums = np.array(sums)
        tot = sums.sum()
        mean_cov = tot / max(1, length)
        ok = mean_cov >= min_mean_cov and tot > 0
        row = {"id": fid, f"{prefix}length": length, f"{prefix}polysome_total": tot,
               f"{prefix}polysome_mean_coverage": mean_cov,
               f"{prefix}MRL": float((np.array(ks) * sums).sum() / tot) if ok else np.nan,
               f"{prefix}frac_polysomal_ge2": float(sums[np.array(ks) >= 2].sum() / tot) if ok else np.nan}
        for k, v in zip(ks, sums):
            row[f"{prefix}frac_k{k}"] = v / tot if tot > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fraction-tracks", required=True, help="TSV: n_ribosomes, strand (+/-), path")
    ap.add_argument("--transcripts-gtf", required=True)
    ap.add_argument("--orfs", required=True)
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--chromosomes", default="all")
    ap.add_argument("--min-mean-coverage", type=float, default=1.0,
                    help="min summed (over fractions) polysome coverage per nt to report an MRL")
    ap.add_argument("--n-cores", type=int, default=4)
    args = ap.parse_args()

    tt = pd.read_csv(args.fraction_tracks, sep="\t")
    ks = sorted(tt["n_ribosomes"].astype(int).unique())
    chroms = set(read_chromosomes(args.chromosomes, args.chrom_sizes))

    ex = read_gtf(args.transcripts_gtf, features={"exon"}, keys=("transcript_id",))
    tx_recs = [(tid, g["chrom"].iat[0], g["strand"].iat[0], list(zip(g["start0"], g["end1"])))
               for tid, g in ex.groupby("transcript_id")]
    orfs = pd.read_csv(args.orfs, sep="\t")
    orf_recs = [(r.orf_id, r.chrom, r.strand, [tuple(map(int, b.split("-"))) for b in r.orf_blocks.split(";")])
                for r in orfs.itertuples()]
    ref = read_gtf(args.reference_gtf, features={"exon"}, keys=("transcript_id", "transcript_type", "tag",
                                                               "gene_name"))
    ref = ref.loc[ref["chrom"].isin(chroms) & (ref["transcript_type"] == "protein_coding")
                  & ref["tag"].fillna("").str.contains("Ensembl_canonical")]
    names = ref.groupby("transcript_id")["gene_name"].first()
    ref_recs = [(tid, g["chrom"].iat[0], g["strand"].iat[0], list(zip(g["start0"], g["end1"])))
                for tid, g in ref.groupby("transcript_id")]

    with mp.Pool(args.n_cores, initializer=init_worker, initargs=(tt,)) as pool:
        tx_res = pool.map(feature_sums, tx_recs, chunksize=50)
        orf_res = pool.map(feature_sums, orf_recs, chunksize=200)
        ref_res = pool.map(feature_sums, ref_recs, chunksize=50)

    tx_df = summarize(tx_res, ks, args.min_mean_coverage).rename(columns={"id": "transcript_id"})
    orf_df = summarize(orf_res, ks, args.min_mean_coverage, prefix="orf_").rename(columns={"id": "orf_id"})
    ref_df = summarize(ref_res, ks, args.min_mean_coverage).rename(columns={"id": "transcript_id"})
    ref_df["gene_name"] = ref_df["transcript_id"].map(names)

    ref_mrl = np.sort(ref_df["MRL"].dropna().values)

    def pct(x):
        return np.searchsorted(ref_mrl, x, side="right") / len(ref_mrl) if len(ref_mrl) and np.isfinite(x) else np.nan

    tx_df["MRL_pct_vs_mRNA"] = tx_df["MRL"].map(pct)
    orf_df["orf_MRL_pct_vs_mRNA"] = orf_df["orf_MRL"].map(pct)
    tx_df.to_csv("polysome.transcripts.tsv", sep="\t", index=False)
    orf_df.to_csv("polysome.orfs.tsv", sep="\t", index=False)
    ref_df.to_csv("polysome.reference_mRNAs.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    bins = np.linspace(1, max(ks), 46)
    for lab, v in [("canonical mRNAs", ref_df["MRL"]), ("lincRNA isoforms", tx_df["MRL"]),
                   ("lincRNA ORFs (ORF blocks)", orf_df["orf_MRL"])]:
        v = v.dropna()
        axes[0].hist(v, bins=bins, density=True, histtype="step", label=f"{lab} (n={len(v)})")
    axes[0].set(xlabel="mean ribosome load", ylabel="density", title="MRL distributions")
    axes[0].legend(fontsize=8)
    for lab, d in [("canonical mRNAs", ref_df), ("lincRNA isoforms", tx_df)]:
        axes[1].scatter(np.log10(d["polysome_mean_coverage"] + 1e-3), d["MRL"], s=3, alpha=0.4, label=lab)
    axes[1].set(xlabel="log10 mean polysome coverage (sum over fractions)", ylabel="MRL")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("polysome_mrl.pdf")
    print(f"mRNA reference: n={len(ref_mrl)}, median MRL={np.median(ref_mrl) if len(ref_mrl) else np.nan:.2f}")
    print(tx_df[["MRL", "MRL_pct_vs_mRNA"]].describe().to_string())


if __name__ == "__main__":
    main()
