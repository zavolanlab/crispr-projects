#!/usr/bin/env python
"""
Genome-wide, strand-specific P-site counts from a genomic Ribo-seq BAM.

Each read of a calibrated class (aligned length x 5' soft-clip, see
riboseq_calibrate_offsets.py) contributes one P-site at the position `offset`
nucleotides downstream of its aligned 5' end (splice-aware). Spliced reads are kept
only if all their introns are in the junction whitelist.

Outputs (sorted bedGraph, one line per position with >= 1 P-site):
  <prefix>.plus.bedGraph, <prefix>.minus.bedGraph                 all reads
  <prefix>.unique.plus.bedGraph, <prefix>.unique.minus.bedGraph   MAPQ == 255 only
  <prefix>.psite_stats.tsv                                         read accounting
"""
import argparse
import collections

import numpy as np
import pandas as pd
import pysam

from lincorfz_utils import read_chromosomes
from ribo_reads import (clip5_of, has_skip, junctions_ok, load_offsets, load_whitelist,
                        psite_genomic)


def write_bedgraph(fh, chrom, arr):
    nz = np.flatnonzero(arr)
    if len(nz):
        pd.DataFrame({"c": chrom, "s": nz, "e": nz + 1, "v": arr[nz]}).to_csv(
            fh, sep="\t", header=False, index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--offsets", required=True)
    ap.add_argument("--whitelist", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--chromosomes", default="all")
    ap.add_argument("--max-clip5", type=int, default=5)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    chroms = read_chromosomes(args.chromosomes, args.chrom_sizes)
    sizes = pd.read_csv(args.chrom_sizes, sep="\t", header=None, index_col=0)[1].to_dict()
    offsets = load_offsets(args.offsets)
    if not offsets:
        raise SystemExit("no calibrated read classes passed QC")
    whitelist = load_whitelist(args.whitelist, chroms)
    stats = collections.Counter()

    outs = {k: open(f"{args.prefix}.{k}.bedGraph", "w")
            for k in ["plus", "minus", "unique.plus", "unique.minus"]}
    bam = pysam.AlignmentFile(args.bam)
    bam_chroms = set(bam.references)
    for chrom in chroms:
        if chrom not in bam_chroms:
            continue
        n = sizes[chrom]
        arrs = {k: np.zeros(n, dtype=np.int32) for k in outs}
        cwl = whitelist.get(chrom, set())
        for r in bam.fetch(chrom):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            stats["reads"] += 1
            off = offsets.get((r.query_alignment_length, clip5_of(r, args.max_clip5)))
            if off is None:
                stats["skipped_uncalibrated_class"] += 1
                continue
            if has_skip(r) and not junctions_ok(r, cwl):
                stats["skipped_unsupported_junction"] += 1
                continue
            p = psite_genomic(r, off)
            if p is None or p < 0 or p >= n:
                stats["skipped_psite_out_of_alignment"] += 1
                continue
            strand = "minus" if r.is_reverse else "plus"
            arrs[strand][p] += 1
            stats["psites"] += 1
            if r.mapping_quality == 255:
                arrs["unique." + strand][p] += 1
                stats["psites_unique"] += 1
        for k, fh in outs.items():
            write_bedgraph(fh, chrom, arrs[k])
        stats[f"psites_{chrom}"] = int(arrs["plus"].sum() + arrs["minus"].sum())
        print(chrom, stats[f"psites_{chrom}"], flush=True)
    bam.close()
    for fh in outs.values():
        fh.close()
    pd.Series(stats, name="value").rename_axis("metric").reset_index().to_csv(
        f"{args.prefix}.psite_stats.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
