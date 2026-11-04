#!/usr/bin/env python
"""
Calibrate Ribo-seq P-site offsets on annotated start codons, per read class
(aligned length x 5' soft-clip), using genomic alignments.

For every Ensembl_canonical protein-coding transcript with a complete CDS:
  * reads whose aligned 5' end falls in [-40, +10] nt (transcript coordinates)
    around the first nucleotide of the start codon give the distance histogram
    H_class[d], d = start - 5'end;
  * reads whose 5' end falls in the CDS body give the reading-frame histogram
    F_class[r], r = (5'end - start) mod 3.

The offset of a class is, among the offsets that place the P-site in the dominant CDS frame,
the one closest to the start-codon peak (the d with the largest H_class[d]; search range
--min-offset..--max-offset).
A class passes when it has enough reads and >= --min-frame0 of its CDS P-sites in frame 0.

Outputs: <prefix>.offsets.tsv, <prefix>.metagene.tsv, <prefix>.calibration.pdf
"""
import argparse
import collections

import matplotlib
import numpy as np
import pandas as pd
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from lincorfz_utils import (TxMapper, blocks_from_positions, read_chromosomes, read_gtf,  # noqa: E402
                            tx_genomic_positions)
from ribo_reads import clip5_of, five_prime_genomic, has_skip, junctions_ok, load_whitelist  # noqa: E402

D_MIN, D_MAX = -10, 40


def canonical_cds_transcripts(gtf_path, chroms, min_cds_len):
    df = read_gtf(gtf_path, features={"exon", "start_codon", "stop_codon"},
                  keys=("transcript_id", "transcript_type", "tag"))
    df = df.loc[df["chrom"].isin(chroms)]
    tags = df.groupby("transcript_id")["tag"].first().fillna("")
    ttype = df.groupby("transcript_id")["transcript_type"].first()
    keep = tags.index[tags.str.contains("Ensembl_canonical")
                      & ~tags.str.contains("cds_start_NF|cds_end_NF|mRNA_start_NF|mRNA_end_NF")
                      & (ttype.reindex(tags.index) == "protein_coding")]
    df = df.loc[df["transcript_id"].isin(set(keep))]
    out = []
    for tid, g in df.groupby("transcript_id"):
        ex = g.loc[g["feature"] == "exon"].sort_values("start0")
        sc = g.loc[g["feature"] == "start_codon"]
        st = g.loc[g["feature"] == "stop_codon"]
        if ex.empty or sc.empty or st.empty:
            continue
        chrom, strand = ex["chrom"].iat[0], ex["strand"].iat[0]
        exons = list(zip(ex["start0"], ex["end1"]))
        m = TxMapper(exons, strand)
        if strand == "+":
            s_g, e_g = sc["start0"].min(), st["end1"].max() - 1
        else:
            s_g, e_g = sc["end1"].max() - 1, st["start0"].min()
        s_tx, e_tx = m.to_tx(s_g), m.to_tx(e_g)
        if s_tx < 0 or e_tx < 0 or (e_tx - s_tx + 1) % 3 != 0 or (e_tx - s_tx + 1) < min_cds_len:
            continue
        out.append((tid, chrom, strand, exons, m, s_tx, e_tx))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--whitelist", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--chromosomes", default="all")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--min-len", type=int, default=25)
    ap.add_argument("--max-len", type=int, default=36)
    ap.add_argument("--max-clip5", type=int, default=5)
    ap.add_argument("--min-offset", type=int, default=8)
    ap.add_argument("--max-offset", type=int, default=20)
    ap.add_argument("--min-start-reads", type=int, default=50)
    ap.add_argument("--min-cds-reads", type=int, default=300)
    ap.add_argument("--min-frame0", type=float, default=0.45)
    ap.add_argument("--cds-window", type=int, default=600, help="CDS body window length (nt)")
    ap.add_argument("--min-cds-len", type=int, default=300)
    args = ap.parse_args()

    chroms = read_chromosomes(args.chromosomes, args.chrom_sizes)
    whitelist = load_whitelist(args.whitelist, chroms)
    txs = canonical_cds_transcripts(args.gtf, chroms, args.min_cds_len)
    print(f"{len(txs)} canonical CDS transcripts used for calibration")

    H = collections.defaultdict(lambda: np.zeros(D_MAX - D_MIN + 1, dtype=np.int64))
    F = collections.defaultdict(lambda: np.zeros(3, dtype=np.int64))
    S = collections.defaultdict(lambda: np.zeros(121, dtype=np.int64))  # 5' end rel. to stop: -90..+30

    bam = pysam.AlignmentFile(args.bam)
    for tid, chrom, strand, exons, m, s_tx, e_tx in txs:
        tx_pos = tx_genomic_positions(exons, strand)
        win_lo, win_hi = max(0, s_tx + D_MIN - 1), min(m.length, s_tx + D_MAX + 1)
        cds_lo = s_tx + 60
        cds_hi = min(s_tx + 60 + args.cds_window, e_tx - 60)
        stop_lo, stop_hi = max(0, e_tx - 2 - 90), min(m.length, e_tx - 2 + 31)
        regions = [(win_lo, win_hi), (stop_lo, stop_hi)]
        if cds_hi > cds_lo:
            regions.append((cds_lo, cds_hi))
        seen = set()
        cwl = whitelist.get(chrom, set())
        for lo, hi in regions:
            for gs, ge in blocks_from_positions(tx_pos[lo:hi]):
                for r in bam.fetch(chrom, gs, ge):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary:
                        continue
                    if r.is_reverse != (strand == "-"):
                        continue
                    key = (r.query_name, r.reference_start)
                    if key in seen:
                        continue
                    seen.add(key)
                    alen = r.query_alignment_length
                    if alen < args.min_len or alen > args.max_len:
                        continue
                    if has_skip(r) and not junctions_ok(r, cwl):
                        continue
                    f_tx = m.to_tx(five_prime_genomic(r))
                    if f_tx < 0:
                        continue
                    cls = (alen, clip5_of(r, args.max_clip5))
                    d = s_tx - f_tx
                    if D_MIN <= d <= D_MAX:
                        H[cls][d - D_MIN] += 1
                    if cds_lo <= f_tx < cds_hi:
                        F[cls][(f_tx - s_tx) % 3] += 1
                    ds = f_tx - (e_tx - 2)
                    if -90 <= ds <= 30:
                        S[cls][ds + 90] += 1
    bam.close()

    rows = []
    offsets_range = np.arange(args.min_offset, args.max_offset + 1)
    for cls in sorted(set(H) | set(F)):
        h, f = H[cls], F[cls]
        n_cds = int(f.sum())
        r_best = int(np.argmax(f)) if n_cds else 0
        allowed = [o for o in offsets_range if (r_best + o) % 3 == 0]
        hv = np.array([h[o - D_MIN] for o in offsets_range])
        o_unconstrained = int(offsets_range[np.argmax(hv)]) if hv.sum() else -1
        # frame-consistent offset closest to the start-codon peak (ties: higher peak)
        o_best = (min(allowed, key=lambda o: (abs(o - o_unconstrained), -h[o - D_MIN]))
                  if allowed and o_unconstrained >= 0 else -1)
        n_start = int(hv.sum())
        frames = [f[(fr - o_best) % 3] for fr in range(3)] if n_cds else [0, 0, 0]
        frame0 = frames[0] / n_cds if n_cds else 0.0
        peak = h[o_best - D_MIN] if o_best >= 0 else 0
        peak_ratio = peak / max(1.0, np.median(hv)) if n_start else 0.0
        ok = (n_start >= args.min_start_reads and n_cds >= args.min_cds_reads
              and frame0 >= args.min_frame0 and frames[0] == max(frames))
        rows.append(dict(aligned_len=cls[0], clip5=cls[1], offset=o_best,
                         offset_unconstrained=o_unconstrained, n_start_window=n_start, n_cds=n_cds,
                         frame0=frames[0], frame1=frames[1], frame2=frames[2],
                         frac_frame0=round(frame0, 4), start_peak_ratio=round(peak_ratio, 3),
                         start_peak_agrees=(o_best == o_unconstrained), **{"pass": ok}))
    res = pd.DataFrame(rows)
    res.to_csv(f"{args.prefix}.offsets.tsv", sep="\t", index=False)

    # metagene of P-sites (passing classes) around start and stop codons
    start_meta = np.zeros(D_MAX - D_MIN + 1 + args.max_offset, dtype=np.int64)
    stop_meta = np.zeros(121 + args.max_offset, dtype=np.int64)
    rel_start = np.arange(-D_MAX, -D_MIN + 1)  # P-site rel. to start for offset 0
    for r in res.loc[res["pass"]].itertuples():
        cls = (r.aligned_len, r.clip5)
        # P-site relative to start: o - d  (d indexes H); write into array shifted by offset
        h = H[cls][::-1]  # index i -> d = D_MAX - i -> relative 5' end = -d
        start_meta[r.offset:r.offset + len(h)] += h
        stop_meta[r.offset:r.offset + 121] += S[cls]
    meta = pd.concat([
        pd.DataFrame({"anchor": "start", "position": np.arange(len(start_meta)) + rel_start[0],
                      "psites": start_meta}),
        pd.DataFrame({"anchor": "stop", "position": np.arange(len(stop_meta)) - 90, "psites": stop_meta}),
    ])
    meta.to_csv(f"{args.prefix}.metagene.tsv", sep="\t", index=False)

    # QC figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))
    top = res.loc[res["pass"]].sort_values("n_cds", ascending=False).head(8)
    for r in top.itertuples():
        h = H[(r.aligned_len, r.clip5)]
        axes[0].plot(np.arange(D_MIN, D_MAX + 1), h / max(1, h.sum()),
                     label=f"len {r.aligned_len}, clip5 {r.clip5}: off {r.offset}")
    axes[0].set(xlabel="start codon - aligned 5' end (nt)", ylabel="fraction of reads",
                title="5' ends around start codons")
    axes[0].legend(fontsize=7)
    for anchor, ax in zip(["start", "stop"], axes[1:]):
        mm = meta.loc[(meta["anchor"] == anchor) & meta["position"].between(-40, 40)]
        cols = ["tab:red" if p % 3 == 0 else ("tab:blue" if p % 3 == 1 else "tab:green")
                for p in mm["position"]]
        ax.bar(mm["position"], mm["psites"], color=cols, width=0.9)
        ax.set(xlabel=f"P-site position relative to {anchor} codon (nt)", ylabel="P-sites",
               title=f"P-site metagene, {anchor}")
    fig.tight_layout()
    fig.savefig(f"{args.prefix}.calibration.pdf")
    print(res.loc[res["pass"]].to_string())


if __name__ == "__main__":
    main()
