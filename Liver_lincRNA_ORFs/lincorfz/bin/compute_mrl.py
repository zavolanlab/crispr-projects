#!/usr/bin/env python
"""
Per-nucleotide mean ribosome load (MRL) from polysome-fraction RNA-seq coverage.

For every genomic position x on one strand of one chromosome:

    c_k(x)  = mean over replicates of the (spike-in normalized) coverage in the
              fraction carrying k ribosomes (k = 1 for 80S, 2 for disome, ... 10 for decasome)
    T(x)    = sum_k c_k(x)                       ("polysomal coverage")
    MRL(x)  = sum_k k * c_k(x) / T(x)            (only where T(x) >= --min-total-coverage)

Replicate means use the replicates available for each fraction, so a fraction that
is missing in one replicate does not bias the MRL. Optionally, an MRL is also
computed for each replicate that has every fraction.

Outputs (bedGraph, chromosome names as in --chrom, i.e. 'chr' notation):
    <chrom>.<strand>.mrl.bedGraph
    <chrom>.<strand>.total.bedGraph
    <chrom>.<strand>.frac<k>.bedGraph           (replicate-mean coverage per fraction)
    <chrom>.<strand>.mrl_rep<r>.bedGraph         (with --per-replicate)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyBigWig

from lincorfz_utils import bigwig_chrom_name


def write_runs(fh, chrom, start, values, decimals, skip_zero=True, valid=None):
    """Write run-length encoded bedGraph lines for values[i] at position start+i."""
    if valid is None:
        valid = np.ones(len(values), dtype=bool)
    vals = np.round(values, decimals)
    if skip_zero:
        valid = valid & (vals != 0)
    if not valid.any():
        return
    # a run breaks where validity or value changes
    change = np.empty(len(vals), dtype=bool)
    change[0] = True
    change[1:] = (vals[1:] != vals[:-1]) | (valid[1:] != valid[:-1])
    run_starts = np.flatnonzero(change)
    run_ends = np.append(run_starts[1:], len(vals))
    keep = valid[run_starts]
    fmt = "%s\t%d\t%d\t%." + str(decimals) + "f\n"
    for s, e, v in zip(run_starts[keep], run_ends[keep], vals[run_starts[keep]]):
        fh.write(fmt % (chrom, start + s, start + e, v))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True,
                    help="TSV with columns: path, fraction, n_ribosomes, replicate, strand (plus/minus)")
    ap.add_argument("--chrom", required=True, help="chromosome in 'chr' notation, e.g. chr20")
    ap.add_argument("--strand", required=True, choices=["plus", "minus"])
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--min-ribosomes", type=int, default=1)
    ap.add_argument("--max-ribosomes", type=int, default=10)
    ap.add_argument("--min-total-coverage", type=float, default=10.0)
    ap.add_argument("--chunk-size", type=int, default=5_000_000)
    ap.add_argument("--per-replicate", action="store_true")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    tab = pd.read_csv(args.table, sep="\t")
    tab = tab.loc[(tab["strand"] == args.strand)
                  & tab["n_ribosomes"].between(args.min_ribosomes, args.max_ribosomes)].copy()
    if tab.empty:
        sys.exit(f"no polysome files for strand {args.strand}")
    ks = sorted(tab["n_ribosomes"].unique())
    chrom_len = int(pd.read_csv(args.chrom_sizes, sep="\t", header=None, index_col=0)[1][args.chrom])

    handles = {}
    for p in tab["path"]:
        handles[p] = pyBigWig.open(p)
    tab["bw_chrom"] = [bigwig_chrom_name(handles[p], args.chrom) for p in tab["path"]]

    reps_complete = []
    if args.per_replicate:
        for r, g in tab.groupby("replicate"):
            if set(g["n_ribosomes"]) == set(ks):
                reps_complete.append(r)
            else:
                print(f"[info] replicate {r} lacks fractions {sorted(set(ks) - set(g['n_ribosomes']))}; "
                      f"no per-replicate MRL", file=sys.stderr)

    os.makedirs(args.outdir, exist_ok=True)
    pre = os.path.join(args.outdir, f"{args.chrom}.{args.strand}")
    out = {"mrl": open(pre + ".mrl.bedGraph", "w"), "total": open(pre + ".total.bedGraph", "w")}
    for k in ks:
        out[f"frac{k}"] = open(pre + f".frac{k}.bedGraph", "w")
    for r in reps_complete:
        out[f"mrl_rep{r}"] = open(pre + f".mrl_rep{r}.bedGraph", "w")

    def fetch(row, s, e):
        if row["bw_chrom"] is None:
            return np.zeros(e - s, dtype=np.float64)
        v = handles[row["path"]].values(row["bw_chrom"], s, e, numpy=True)
        return np.nan_to_num(v, nan=0.0).astype(np.float64)

    for s in range(0, chrom_len, args.chunk_size):
        e = min(s + args.chunk_size, chrom_len)
        cov = {}  # (k, replicate) -> array
        for _, row in tab.iterrows():
            cov[(row["n_ribosomes"], row["replicate"])] = fetch(row, s, e)
        mean_k = {}
        for k in ks:
            arrs = [v for (kk, _), v in cov.items() if kk == k]
            mean_k[k] = np.mean(arrs, axis=0)
        total = np.sum([mean_k[k] for k in ks], axis=0)
        weighted = np.sum([k * mean_k[k] for k in ks], axis=0)
        valid = total >= args.min_total_coverage
        mrl = np.zeros_like(total)
        mrl[valid] = weighted[valid] / total[valid]
        write_runs(out["mrl"], args.chrom, s, mrl, 2, skip_zero=False, valid=valid)
        write_runs(out["total"], args.chrom, s, total, 2)
        for k in ks:
            write_runs(out[f"frac{k}"], args.chrom, s, mean_k[k], 2)
        for r in reps_complete:
            tot_r = np.sum([cov[(k, r)] for k in ks], axis=0)
            w_r = np.sum([k * cov[(k, r)] for k in ks], axis=0)
            val_r = tot_r >= args.min_total_coverage
            mrl_r = np.zeros_like(tot_r)
            mrl_r[val_r] = w_r[val_r] / tot_r[val_r]
            write_runs(out[f"mrl_rep{r}"], args.chrom, s, mrl_r, 2, skip_zero=False, valid=val_r)

    for fh in out.values():
        fh.close()
    for h in handles.values():
        h.close()


if __name__ == "__main__":
    main()
