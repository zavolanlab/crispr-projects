#!/usr/bin/env python
"""
Evolutionary conservation of ORFs, following the quantification used in the HCC Ribo-seq
notebook (HCC_RiboSeq_paper.ipynb, "Specify association jobs accounting for track specifics"):

  1. every genomic segment of an ORF gets the frame-specific PhyloCSF track matching its
     codon phase: with cum_len = ORF nucleotides upstream of the segment (transcript order)
     and offset = (3 - cum_len % 3) % 3,
        '+' strand: frame = ((start0 + offset) % 3) + 1                 -> PhyloCSF+<frame>
        '-' strand: frame = ((chrom_size - (end1 - offset)) % 3) + 1    -> PhyloCSF-<frame>
  2. segments overlapping annotated protein-coding CDSs are removed with
     `bedtools subtract` (strand-unaware, as in the notebook), so that conservation of
     overlapping coding genes does not leak into the ORF scores;
  3. per ORF, mean and median over all remaining positions with defined values.

The same subtracted segments are scored with PhyloCSF_Power, PhyloP100way and PhyloP30way.

Output: conservation.orfs.tsv
"""
import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd
import pyBigWig

from lincorfz_utils import bigwig_chrom_name, read_gtf


def phylocsf_track(strand, start0, end1, cum_len, chrom_size):
    offset = (3 - (cum_len % 3)) % 3
    if strand == "+":
        return f"PhyloCSF+{((start0 + offset) % 3) + 1}"
    return f"PhyloCSF-{((chrom_size - (end1 - offset)) % 3) + 1}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orfs", required=True)
    ap.add_argument("--reference-gtf", required=True, help="annotation providing CDS features to subtract")
    ap.add_argument("--conservation-dir", required=True,
                    help="directory with PhyloCSF+1..3/-1..3, PhyloCSF_Power, PhyloP100way, PhyloP30way .bw")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--extra-tracks", default="PhyloCSF_Power,PhyloP100way,PhyloP30way")
    args = ap.parse_args()

    sizes = pd.read_csv(args.chrom_sizes, sep="\t", header=None, index_col=0)[1].to_dict()
    orfs = pd.read_csv(args.orfs, sep="\t")
    segs = []
    for r in orfs.itertuples():
        blocks = [tuple(map(int, b.split("-"))) for b in r.orf_blocks.split(";")]
        order = blocks if r.strand == "+" else blocks[::-1]
        cum = 0
        for s, e in order:
            segs.append((r.chrom, s, e, r.orf_id, phylocsf_track(r.strand, s, e, cum, sizes[r.chrom])))
            cum += e - s
    segs = pd.DataFrame(segs, columns=["chrom", "start0", "end1", "orf_id", "track"]).sort_values(
        ["chrom", "start0"])

    cds = read_gtf(args.reference_gtf, features={"CDS"}, keys=("transcript_id",))
    tmp = tempfile.mkdtemp(dir=".")
    a, b = os.path.join(tmp, "orf_segments.bed"), os.path.join(tmp, "cds.bed")
    segs.to_csv(a, sep="\t", header=False, index=False)
    cds[["chrom", "start0", "end1"]].sort_values(["chrom", "start0"]).to_csv(b, sep="\t", header=False,
                                                                              index=False)
    out = subprocess.run(["bedtools", "subtract", "-a", a, "-b", b], check=True, capture_output=True,
                         text=True).stdout
    kept = pd.DataFrame([l.split("\t") for l in out.splitlines()],
                        columns=["chrom", "start0", "end1", "orf_id", "track"])
    kept[["start0", "end1"]] = kept[["start0", "end1"]].astype(int)
    subprocess.run(["rm", "-rf", tmp], check=True)

    handles = {}

    def bw(name):
        if name not in handles:
            p = os.path.join(args.conservation_dir, f"{name}.bw")
            handles[name] = pyBigWig.open(p) if os.path.exists(p) else None
        return handles[name]

    track_sets = {"PhyloCSF": None}
    for t in [x for x in args.extra_tracks.split(",") if x]:
        track_sets[t] = t
    vals = {k: {} for k in track_sets}
    for r in kept.itertuples():
        for key, tname in track_sets.items():
            h = bw(r.track if tname is None else tname)
            if h is None:
                continue
            c = bigwig_chrom_name(h, r.chrom)
            if c is None:
                continue
            try:
                v = np.asarray(h.values(c, r.start0, r.end1), dtype=float)
            except RuntimeError:
                continue
            vals[key].setdefault(r.orf_id, []).append(v[np.isfinite(v)])
    res = pd.DataFrame({"orf_id": orfs["orf_id"]})
    covered = kept.assign(len=kept["end1"] - kept["start0"]).groupby("orf_id")["len"].sum()
    res["cons_scored_nt"] = res["orf_id"].map(covered).fillna(0).astype(int)
    for key in track_sets:
        d = {k: np.concatenate(v) for k, v in vals[key].items()}
        res[f"{key}_Median"] = res["orf_id"].map({k: np.median(v) if len(v) else np.nan for k, v in d.items()})
        res[f"{key}_Mean"] = res["orf_id"].map({k: np.mean(v) if len(v) else np.nan for k, v in d.items()})
    res.to_csv("conservation.orfs.tsv", sep="\t", index=False)
    for h in handles.values():
        if h is not None:
            h.close()
    print(res.describe().T.to_string())


if __name__ == "__main__":
    main()
