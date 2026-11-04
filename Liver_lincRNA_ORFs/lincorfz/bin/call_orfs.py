#!/usr/bin/env python
"""
Call translated ORFs in de novo lincRNA isoforms from genomic Ribo-seq P-sites.

For every selected isoform, all ORFs (start codon .. first in-frame stop, >= --min-aa codons)
are enumerated on the spliced sequence. ORFs sharing a stop codon are reduced to one start:
the most upstream start whose segment up to the next candidate start carries in-frame P-site
support (in-frame P-sites >= 1 and > each off-frame); otherwise the most upstream ATG.

Per ORF and per P-site track (each replicate, pooled, pooled unique reads):
  * P-sites in frames 0/1/2 over sense codons, in-frame fraction, codons with in-frame P-sites;
  * periodicity p-value as in the HCC Ribo-seq paper: one-sided Wilcoxon rank-sum tests of the
    per-codon in-frame P-sites against each off-frame, combined with Fisher's method;
for the pooled track additionally
  * island score: log2 of the P-site density in the ORF over the density in the rest of the isoform;
  * fraction of the isoform's P-sites inside the ORF, fraction from uniquely mapped reads;
  * start codon, Kozak context (-3 purine, +4 G), amino-acid sequence, genomic CDS blocks.

ORFs with identical genomic CDS structure found in several isoforms are collapsed; the
representative isoform is the most expressed one (TPM_max).

With --length-tracks, P-sites in frames 0/1/2 of every ORF are also counted per sample and aligned
read length (input of the Bayesian frame test, bayes_frame_test.py).

Outputs: orfs.tsv, orf_isoforms.tsv, transcript_ribo.tsv, orf_length_frame_counts.tsv
"""
import argparse
import multiprocessing as mp

import numpy as np
import pandas as pd
import pyBigWig
import pysam
from scipy import stats

from lincorfz_utils import (STOP_CODONS, bigwig_along_tx, bigwig_chrom_name, blocks_from_positions,
                            read_gtf, spliced_sequence, translate, tx_genomic_positions)

G = {}  # per-worker globals


def init_worker(genome, track_table, min_psites, length_table):
    G["fasta"] = pysam.FastaFile(genome)
    G["tracks"] = {}
    for r in track_table.itertuples():
        G["tracks"][(r.label, r.strand)] = pyBigWig.open(r.path)
    G["labels"] = list(dict.fromkeys(track_table["label"]))
    G["min_psites"] = min_psites
    # per-read-length raw P-site tracks: (sample, length, strand) -> bigwig
    G["ltracks"] = {}
    if length_table is not None:
        for r in length_table.itertuples():
            G["ltracks"][(r.sample, int(r.length), r.strand)] = pyBigWig.open(r.path)


def periodicity_p(codon_mat):
    """Fisher-combined one-sided Wilcoxon rank-sum p (frame 0 > frame 1, frame 0 > frame 2)."""
    f0, f1, f2 = codon_mat[:, 0], codon_mat[:, 1], codon_mat[:, 2]
    if f0.sum() == 0 or len(f0) < 2:
        return 1.0
    ps = []
    for other in (f1, f2):
        try:
            p = stats.mannwhitneyu(f0, other, alternative="greater", method="asymptotic").pvalue
        except ValueError:
            p = 1.0
        ps.append(1.0 if not np.isfinite(p) else max(p, 1e-300))
    chi = -2.0 * np.sum(np.log(ps))
    return float(stats.chi2.sf(chi, 2 * len(ps)))


def choose_start(starts, stop_i, inframe_ok):
    for k, s in enumerate(starts):
        nxt = starts[k + 1] if k + 1 < len(starts) else stop_i
        if inframe_ok(s, nxt):
            return s
    return None


def process_transcript(rec):
    tid, chrom, strand, exons, min_aa, start_set, near_cognate = rec
    fasta, tracks, labels = G["fasta"], G["tracks"], G["labels"]
    seq = spliced_sequence(fasta, chrom, exons, strand)
    L = len(seq)
    gpos = tx_genomic_positions(exons, strand)
    ps = {}
    for lab in labels:
        bw = tracks[(lab, strand)]
        ps[lab] = bigwig_along_tx(bw, bigwig_chrom_name(bw, chrom), exons, strand)
    pooled = ps["pooled"]
    tx_total = float(pooled.sum())
    tx_row = dict(transcript_id=tid, tx_length=L, tx_psites=tx_total,
                  tx_psites_unique=float(ps["pooled_unique"].sum()))
    out = []
    if tx_total == 0:
        return out, tx_row, []

    def frame_sums(a, s, e):
        seg = a[s:e]
        n = (e - s) // 3
        m = seg[: n * 3].reshape(n, 3)
        return m

    def inframe_ok(s, e):
        if e - s < 3:
            return False
        m = frame_sums(pooled, s, e).sum(axis=0)
        return m[0] >= 1 and m[0] > m[1] and m[0] > m[2]

    all_starts = start_set | near_cognate
    for frame in range(3):
        cand = []
        for i in range(frame, L - 2, 3):
            cod = seq[i:i + 3]
            if cod in STOP_CODONS:
                if cand:
                    atg = [s for s in cand if seq[s:s + 3] in start_set]
                    supported_nc = [s for s in cand if seq[s:s + 3] in near_cognate
                                    and inframe_ok(s, s + 30 if s + 30 < i else i)]
                    pool = sorted(set(atg) | set(supported_nc))
                    if pool and (i - pool[0]) // 3 >= min_aa:
                        s = choose_start(pool, i, inframe_ok)
                        if s is None:
                            s = atg[0] if atg else None
                        if s is not None and (i - s) // 3 >= min_aa:
                            out.append((s, i, len(pool), pool[0]))
                cand = []
            elif cod in all_starts:
                cand.append(i)

    rows = []
    for s, stop_i, n_alt, longest_start in out:
        e = stop_i  # sense codons [s, stop_i)
        orf_pooled = float(pooled[s:e].sum())
        if orf_pooled < G["min_psites"]:
            continue
        row = dict(transcript_id=tid, chrom=chrom, strand=strand, tx_start=s, tx_stop_end=stop_i + 3,
                   n_codons=(e - s) // 3, n_alt_starts=n_alt, longest_start=longest_start)
        for lab in labels:
            m = frame_sums(ps[lab], s, e)
            fr = m.sum(axis=0)
            tot = float(fr.sum())
            row[f"{lab}_psites"] = tot
            row[f"{lab}_f0"], row[f"{lab}_f1"], row[f"{lab}_f2"] = map(float, fr)
            row[f"{lab}_frac_inframe"] = fr[0] / tot if tot > 0 else np.nan
            row[f"{lab}_codons_inframe"] = int((m[:, 0] > 0).sum())
            row[f"{lab}_p_periodicity"] = periodicity_p(m) if tot > 0 else 1.0
        out_len = L - (e - s)
        out_ps = tx_total - orf_pooled
        row["island_log2ratio"] = float(np.log2((orf_pooled + 0.5) / (e - s))
                                        - np.log2((out_ps + 0.5) / max(out_len, 30)))
        row["frac_tx_psites_in_orf"] = orf_pooled / tx_total
        row["frac_psites_unique"] = (float(ps["pooled_unique"][s:e].sum()) / orf_pooled) if orf_pooled else np.nan
        row["start_codon"] = seq[s:s + 3]
        row["stop_codon"] = seq[stop_i:stop_i + 3]
        row["kozak_context"] = seq[max(0, s - 6):s + 4]
        row["kozak_score"] = int(s >= 3 and seq[s - 3] in "AG") + int(seq[s + 3:s + 4] == "G")
        row["aa_seq"] = translate(seq[s:stop_i])
        blocks = blocks_from_positions(gpos[s:stop_i + 3])
        row["orf_start0"] = blocks[0][0]
        row["orf_end1"] = blocks[-1][1]
        row["orf_blocks"] = ";".join(f"{a}-{b}" for a, b in blocks)
        rows.append(row)
    # frame counts of every reported ORF per sample and aligned read length
    lrows = []
    if rows and G["ltracks"]:
        for (sample, rl, st), bw in G["ltracks"].items():
            if st != strand:
                continue
            v = bigwig_along_tx(bw, bigwig_chrom_name(bw, chrom), exons, strand)
            if v.sum() == 0:
                continue
            vp = np.minimum(v, 1.0)   # covered positions
            for row in rows:
                s0, e0 = row["tx_start"], row["tx_stop_end"] - 3
                fr = frame_sums(v, s0, e0).sum(axis=0)
                if fr.sum() > 0:
                    fp = frame_sums(vp, s0, e0).sum(axis=0)
                    lrows.append((tid, s0, row["tx_stop_end"], sample, rl, *map(float, fr), *map(float, fp)))
    return rows, tx_row, lrows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--transcripts-gtf", required=True)
    ap.add_argument("--transcripts-tsv", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--psite-tracks", required=True,
                    help="TSV: label, strand (+/-), path; labels must include 'pooled' and 'pooled_unique'")
    ap.add_argument("--length-tracks", default="",
                    help="TSV: sample, length, strand (+/-), path of per-read-length raw P-site bigwigs")
    ap.add_argument("--min-aa", type=int, default=8)
    ap.add_argument("--start-codons", default="ATG")
    ap.add_argument("--near-cognate-codons", default="", help="e.g. CTG,GTG,TTG (empty = ATG only)")
    ap.add_argument("--min-psites", type=float, default=1.0, help="min pooled P-sites to report an ORF")
    ap.add_argument("--n-cores", type=int, default=4)
    args = ap.parse_args()

    tracks = pd.read_csv(args.psite_tracks, sep="\t")
    ex = read_gtf(args.transcripts_gtf, features={"exon"}, keys=("transcript_id",))
    info = pd.read_csv(args.transcripts_tsv, sep="\t").set_index("transcript_id")
    start_set = set(c.strip() for c in args.start_codons.split(",") if c.strip())
    near = set(c.strip() for c in args.near_cognate_codons.split(",") if c.strip()) - start_set
    recs = []
    for tid, g in ex.groupby("transcript_id"):
        g = g.sort_values("start0")
        recs.append((tid, g["chrom"].iat[0], g["strand"].iat[0], list(zip(g["start0"], g["end1"])),
                     args.min_aa, start_set, near))
    print(f"{len(recs)} isoforms")

    ltab = pd.read_csv(args.length_tracks, sep="\t") if args.length_tracks else None
    with mp.Pool(args.n_cores, initializer=init_worker,
                 initargs=(args.genome, tracks, args.min_psites, ltab)) as pool:
        res = pool.map(process_transcript, recs, chunksize=20)
    orfs = pd.DataFrame([r for rows, _, _ in res for r in rows])
    txr = pd.DataFrame([t for _, t, _ in res])
    lcnt = pd.DataFrame([r for _, _, lr in res for r in lr],
                        columns=["transcript_id", "tx_start", "tx_stop_end", "sample", "length", "f0", "f1", "f2",
                                 "pos0", "pos1", "pos2"])
    txr.to_csv("transcript_ribo.tsv", sep="\t", index=False)
    if orfs.empty:
        raise SystemExit("no ORFs with P-site support found")

    orfs["TPM_max"] = orfs["transcript_id"].map(info["TPM_max"])
    orfs["locus_id"] = orfs["transcript_id"].map(info["gene_id"])
    orfs["locus_name"] = orfs["transcript_id"].map(info["locus_name"])
    orfs["orf_key"] = orfs["chrom"] + ":" + orfs["strand"] + ":" + orfs["orf_blocks"]
    iso = orfs[["orf_key", "transcript_id", "TPM_max"]].copy()
    rep = orfs.sort_values(["orf_key", "TPM_max"], ascending=[True, False]).drop_duplicates("orf_key")
    isolist = iso.groupby("orf_key")["transcript_id"].agg(lambda x: ",".join(sorted(set(x))))
    rep = rep.copy()
    rep["isoforms"] = rep["orf_key"].map(isolist)
    rep["n_isoforms"] = rep["isoforms"].str.count(",") + 1
    rep["n_aa"] = rep["aa_seq"].str.len()
    rep = rep.sort_values(["chrom", "orf_start0"]).reset_index(drop=True)
    rep.insert(0, "orf_id", [f"{t}|{s}-{e}" for t, s, e in zip(rep["transcript_id"], rep["tx_start"],
                                                                   rep["tx_stop_end"])])
    rep.to_csv("orfs.tsv", sep="\t", index=False)
    # per-length frame counts of the representative isoform of each ORF (counts are genomic,
    # hence identical for all isoforms sharing the ORF)
    lcnt = lcnt.merge(rep[["transcript_id", "tx_start", "tx_stop_end", "orf_id"]],
                      on=["transcript_id", "tx_start", "tx_stop_end"], how="inner")
    lcnt[["orf_id", "sample", "length", "f0", "f1", "f2", "pos0", "pos1", "pos2"]].to_csv(
        "orf_length_frame_counts.tsv", sep="\t", index=False)
    iso.merge(rep[["orf_key", "orf_id"]], on="orf_key").to_csv("orf_isoforms.tsv", sep="\t", index=False)
    print(f"{len(rep)} unique ORFs with >= {args.min_psites} pooled P-sites")


if __name__ == "__main__":
    main()
