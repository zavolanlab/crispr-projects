#!/usr/bin/env python
"""
Read-length-specific priors on the reading-frame distribution of Ribo-seq P-sites.

For every read length L and dataset (each Ribo-seq sample and the pooled samples), Ensembl_canonical
protein-coding transcripts are cut into ORF-sized windows of --window-codons codons:
  * translated model: windows of the CDS body (from codon 16 to 5 codons before the stop codon;
    initiation and termination peaks excluded), frames relative to the annotated reading frame;
  * untranslated model: windows of the 3'UTR, starting 30 nt after the stop codon, frames relative
    to the (continued) CDS frame -- reads there are not produced by translation of the CDS and
    approximate the frame distribution of non-translational background.
A Dirichlet-multinomial distribution is fitted by maximum likelihood to the per-window P-site count
vectors (windows with >= --min-window-psites P-sites of that length), giving alpha_cds(L) and
alpha_null(L). Fitting both models on windows of the size of a typical lincRNA ORF makes their
dispersion comparable to that of ORF-level counts: reads stacking at a few positions, frequent in
short regions, then inflate the variance of both models in the same way (whole-gene CDS counts
would give an unrealistically tight translated model, and sparse 3'UTR counts a very loose null).

The mean a0 / sum(a) is the expected in-frame fraction of a translated ORF for reads of length L;
the concentration sum(a) captures how much this fraction varies between windows. Lengths with poor
frame preservation get a0 ~ a1 ~ a2 and contribute little evidence in the Bayesian frame test
(bayes_frame_test.py). With --position-cap 1, covered positions are counted instead of reads (the
ORF counts must then be positions too).

Outputs: frame_priors.tsv, frame_counts_by_length.tsv, frame_priors.pdf
"""
import argparse
import multiprocessing as mp

import matplotlib
import numpy as np
import pandas as pd
import pyBigWig
from scipy.optimize import minimize
from scipy.special import gammaln

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from lincorfz_utils import bigwig_along_tx, bigwig_chrom_name, read_chromosomes  # noqa: E402
from riboseq_calibrate_offsets import canonical_cds_transcripts  # noqa: E402

G = {}


def init_worker(tracks, cap, w):
    G["bw"] = {(r.sample, int(r.length), r.strand): pyBigWig.open(r.path) for r in tracks.itertuples()}
    G["cap"] = cap
    G["w"] = w


def windows(v, start, stop, w):
    """(n_windows, 3) frame counts of consecutive windows of w codons in v[start:stop] (start in frame 0)."""
    n = (stop - start) // (3 * w)
    if n <= 0:
        return np.zeros((0, 3))
    return v[start:start + n * 3 * w].reshape(n, w, 3).sum(axis=1)


def count_transcript(rec):
    tid, chrom, strand, exons, s_tx, e_tx = rec
    w = G["w"]
    per_len = {}
    for (sample, L, st), bw in G["bw"].items():
        if st == strand:
            v = bigwig_along_tx(bw, bigwig_chrom_name(bw, chrom), exons, strand)
            per_len.setdefault(L, {})[sample] = v
    out = []
    cds_start, cds_stop = s_tx + 45, e_tx - 2 - 15
    utr_start = s_tx + 3 * int(np.ceil((e_tx + 31 - s_tx) / 3))
    for L, d in per_len.items():
        if G["cap"] > 0:
            d = {s: np.minimum(v, G["cap"]) for s, v in d.items()}
        d = dict(d, pooled=np.sum(list(d.values()), axis=0))
        for ds, v in d.items():
            cw = windows(v, cds_start, cds_stop, w)
            uw = windows(v, utr_start, len(v), w) if utr_start < len(v) else np.zeros((0, 3))
            if cw.sum() + uw.sum() > 0:
                out.append((ds, L, cw[cw.sum(axis=1) > 0], uw[uw.sum(axis=1) > 0]))
    return out


def dm_negll(log_a, X):
    a = np.exp(log_a)
    A = a.sum()
    N = X.sum(axis=1)
    ll = gammaln(A) - gammaln(N + A) + (gammaln(X + a) - gammaln(a)).sum(axis=1)
    return -ll.sum()


def fit_dm(X, max_rows, rng):
    """ML fit of a Dirichlet-multinomial to count vectors X (rows x 3)."""
    if len(X) > max_rows:
        X = X[rng.choice(len(X), max_rows, replace=False)]
    p = (X.sum(axis=0) + 0.5) / (X.sum() + 1.5)
    res = minimize(dm_negll, np.log(p * 10), args=(X,), method="L-BFGS-B", bounds=[(-8, 12)] * 3)
    return np.exp(res.x), res.success


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--length-tracks", required=True, help="TSV: sample, length, strand (+/-), path (raw counts)")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--chromosomes", default="all")
    ap.add_argument("--window-codons", type=int, default=20)
    ap.add_argument("--min-window-psites", type=int, default=5)
    ap.add_argument("--position-cap", type=float, default=0.0,
                    help="cap of P-site counts per position and read length (1 = count covered positions; "
                         "0 = count reads, default)")
    ap.add_argument("--min-windows", type=int, default=50)
    ap.add_argument("--max-windows", type=int, default=30000, help="random subsample for the fit")
    ap.add_argument("--n-cores", type=int, default=4)
    args = ap.parse_args()

    rng = np.random.default_rng(1)
    tracks = pd.read_csv(args.length_tracks, sep="\t")
    chroms = read_chromosomes(args.chromosomes, args.chrom_sizes)
    txs = canonical_cds_transcripts(args.gtf, chroms, 300)
    recs = [(tid, chrom, strand, exons, s_tx, e_tx) for tid, chrom, strand, exons, _m, s_tx, e_tx in txs]
    print(f"{len(recs)} canonical CDS transcripts")
    with mp.Pool(args.n_cores, initializer=init_worker,
                 initargs=(tracks, args.position_cap, args.window_codons)) as pool:
        res = pool.map(count_transcript, recs, chunksize=50)

    cds, utr = {}, {}
    for rows in res:
        for ds, L, cw, uw in rows:
            cds.setdefault((ds, L), []).append(cw)
            utr.setdefault((ds, L), []).append(uw)

    rows, null_conc = [], {}
    keys = sorted(set(cds) | set(utr), key=lambda k: (k[0], k[1]))
    for ds, L in keys:
        C = np.concatenate(cds.get((ds, L), [np.zeros((0, 3))]))
        U = np.concatenate(utr.get((ds, L), [np.zeros((0, 3))]))
        Cw = C[C.sum(axis=1) >= args.min_window_psites]
        Uw = U[U.sum(axis=1) >= args.min_window_psites]
        tc, tu = C.sum(axis=0), U.sum(axis=0)
        row = dict(dataset=ds, length=int(L), window_codons=args.window_codons,
                   cds_counts=float(tc.sum()), utr_counts=float(tu.sum()),
                   cds_frac_f0=tc[0] / max(1, tc.sum()), cds_frac_f1=tc[1] / max(1, tc.sum()),
                   cds_frac_f2=tc[2] / max(1, tc.sum()), utr_frac_f0=tu[0] / max(1, tu.sum()),
                   n_windows_cds=len(Cw), n_windows_utr=len(Uw))
        if len(Cw) >= args.min_windows:
            a, ok = fit_dm(Cw, args.max_windows, rng)
            row.update(alpha_cds_0=a[0], alpha_cds_1=a[1], alpha_cds_2=a[2], fit_cds_ok=ok)
        if len(Uw) >= args.min_windows:
            a, ok = fit_dm(Uw, args.max_windows, rng)
            row.update(alpha_null_0=a[0], alpha_null_1=a[1], alpha_null_2=a[2], fit_null_ok=ok, null_fitted=True)
            null_conc.setdefault(ds, []).append(a.sum())
        rows.append(row)
    pri = pd.DataFrame(rows)
    for c in ["alpha_cds_0", "alpha_cds_1", "alpha_cds_2", "alpha_null_0", "alpha_null_1", "alpha_null_2",
              "null_fitted"]:
        if c not in pri:
            pri[c] = np.nan
    # lengths with too few 3'UTR windows: symmetric null with the dataset's median fitted concentration
    for i, r in pri.iterrows():
        if not r.get("null_fitted", False) or pd.isna(r.get("null_fitted")):
            conc = float(np.median(null_conc[r["dataset"]])) if null_conc.get(r["dataset"]) else 10.0
            pri.loc[i, ["alpha_null_0", "alpha_null_1", "alpha_null_2"]] = conc / 3
            pri.loc[i, "null_fitted"] = False
    a = pri[["alpha_cds_0", "alpha_cds_1", "alpha_cds_2"]]
    pri["prior_inframe_mean"] = a["alpha_cds_0"] / a.sum(axis=1)
    pri["prior_concentration"] = a.sum(axis=1)
    pri["null_concentration"] = pri[["alpha_null_0", "alpha_null_1", "alpha_null_2"]].sum(axis=1)
    pri = pri.sort_values(["dataset", "length"])
    pri.to_csv("frame_priors.tsv", sep="\t", index=False)
    pri[["dataset", "length", "cds_counts", "cds_frac_f0", "cds_frac_f1", "cds_frac_f2", "utr_counts",
         "utr_frac_f0"]].to_csv("frame_counts_by_length.tsv", sep="\t", index=False)

    datasets = list(pri["dataset"].unique())
    fig, axes = plt.subplots(2, len(datasets), figsize=(4.2 * len(datasets), 6.5), squeeze=False, sharex=True)
    unit = "positions" if args.position_cap > 0 else "P-sites"
    for j, d in enumerate(datasets):
        p = pri.loc[pri["dataset"] == d]
        axes[0, j].bar(p["length"], p["cds_counts"] / 1e6, color="grey")
        axes[0, j].set(title=d, ylabel=f"CDS {unit} (millions)")
        for k, col in enumerate(["cds_frac_f0", "cds_frac_f1", "cds_frac_f2"]):
            axes[1, j].plot(p["length"], p[col], marker="o", color=["tab:red", "tab:blue", "tab:green"][k],
                            label=f"CDS frame {k}")
        axes[1, j].plot(p["length"], p["utr_frac_f0"], ls="--", color="k", label="3'UTR frame 0")
        axes[1, j].plot(p["length"], p["prior_inframe_mean"], ls=":", color="tab:red", label="prior mean frame 0")
        axes[1, j].axhline(1 / 3, color="grey", lw=0.5)
        axes[1, j].set(xlabel="aligned read length (nt)", ylabel=f"fraction of {unit}", ylim=(0, 1))
    axes[1, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig("frame_priors.pdf")
    print(pri[["dataset", "length", "cds_counts", "cds_frac_f0", "prior_inframe_mean", "prior_concentration",
               "null_concentration", "n_windows_cds", "n_windows_utr"]].to_string())


if __name__ == "__main__":
    main()
