#!/usr/bin/env python
"""
Genome-wide, strand-specific P-site counts from a genomic Ribo-seq BAM.

Each read contributes one P-site at the position `offset` nucleotides downstream of its
aligned 5' end (splice-aware); offsets are calibrated per read class (aligned length x 5'
soft-clip, see riboseq_calibrate_offsets.py). Spliced reads are kept only if all their introns
are in the junction whitelist.

Two sets of tracks are written:
  * main tracks: reads of classes that passed the calibration (frame-quality) criteria;
  * per-read-length tracks: reads of every class with a defined offset (also classes that failed
    the frame criteria), split by aligned read length. The frame quality of each length is
    assessed downstream (frame priors per length), so that poor lengths can be down-weighted
    instead of being discarded.

Outputs (sorted bedGraph, one line per position with >= 1 P-site):
  <prefix>.plus.bedGraph, <prefix>.minus.bedGraph                 main tracks, all reads
  <prefix>.unique.plus.bedGraph, <prefix>.unique.minus.bedGraph   main tracks, MAPQ == 255 only
  <prefix>.len<L>.plus.bedGraph, <prefix>.len<L>.minus.bedGraph   per aligned read length L
  <prefix>.psite_stats.tsv                                         read accounting
"""
import argparse
import collections
from array import array

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


def write_positions(fh, chrom, positions):
    if len(positions):
        pos, cnt = np.unique(np.frombuffer(positions, dtype=np.int64), return_counts=True)
        pd.DataFrame({"c": chrom, "s": pos, "e": pos + 1, "v": cnt}).to_csv(
            fh, sep="\t", header=False, index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--offsets", required=True)
    ap.add_argument("--whitelist", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--chromosomes", default="all")
    ap.add_argument("--max-clip5", type=int, default=5)
    ap.add_argument("--min-len", type=int, default=25, help="shortest aligned length with its own track")
    ap.add_argument("--max-len", type=int, default=36, help="longest aligned length with its own track")
    ap.add_argument("--min-cds-reads-length-tracks", type=int, default=100,
                    help="classes failing calibration enter the per-length tracks if they have a defined "
                         "offset and at least this many CDS reads")
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    chroms = read_chromosomes(args.chromosomes, args.chrom_sizes)
    sizes = pd.read_csv(args.chrom_sizes, sep="\t", header=None, index_col=0)[1].to_dict()
    offsets = load_offsets(args.offsets)
    offsets_any = load_offsets(args.offsets, include_failed=True,
                               min_cds_reads=args.min_cds_reads_length_tracks)
    if not offsets:
        raise SystemExit("no calibrated read classes passed QC")
    whitelist = load_whitelist(args.whitelist, chroms)
    lengths = list(range(args.min_len, args.max_len + 1))
    stats = collections.Counter()

    outs = {k: open(f"{args.prefix}.{k}.bedGraph", "w")
            for k in ["plus", "minus", "unique.plus", "unique.minus"]}
    louts = {(L, s): open(f"{args.prefix}.len{L}.{s}.bedGraph", "w") for L in lengths for s in ["plus", "minus"]}
    bam = pysam.AlignmentFile(args.bam)
    bam_chroms = set(bam.references)
    for chrom in chroms:
        if chrom not in bam_chroms:
            continue
        n = sizes[chrom]
        arrs = {k: np.zeros(n, dtype=np.int32) for k in outs}
        lpos = {k: array("q") for k in louts}
        cwl = whitelist.get(chrom, set())
        for r in bam.fetch(chrom):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            stats["reads"] += 1
            alen = r.query_alignment_length
            cls = (alen, clip5_of(r, args.max_clip5))
            off_any = offsets_any.get(cls)
            if off_any is None:
                stats["skipped_uncalibrated_class"] += 1
                continue
            if has_skip(r) and not junctions_ok(r, cwl):
                stats["skipped_unsupported_junction"] += 1
                continue
            p = psite_genomic(r, off_any)
            if p is None or p < 0 or p >= n:
                stats["skipped_psite_out_of_alignment"] += 1
                continue
            strand = "minus" if r.is_reverse else "plus"
            if args.min_len <= alen <= args.max_len:
                lpos[(alen, strand)].append(p)
                stats["psites_all_lengths"] += 1
                stats[f"psites_len{alen}"] += 1
            if cls not in offsets:   # failed class: per-length tracks only
                stats["psites_failed_class_length_tracks_only"] += 1
                continue
            arrs[strand][p] += 1
            stats["psites"] += 1
            if r.mapping_quality == 255:
                arrs["unique." + strand][p] += 1
                stats["psites_unique"] += 1
        for k, fh in outs.items():
            write_bedgraph(fh, chrom, arrs[k])
        for k, fh in louts.items():
            write_positions(fh, chrom, lpos[k])
        stats[f"psites_{chrom}"] = int(arrs["plus"].sum() + arrs["minus"].sum())
        print(chrom, stats[f"psites_{chrom}"], flush=True)
    bam.close()
    for fh in list(outs.values()) + list(louts.values()):
        fh.close()
    pd.Series(stats, name="value").rename_axis("metric").reset_index().to_csv(
        f"{args.prefix}.psite_stats.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
