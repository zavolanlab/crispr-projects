#!/usr/bin/env python
"""
Turn peptides identified by FragPipe into the RNA sequences that encode them, as single-end reads
that can be aligned to the genome like Ribo-seq reads.

For every peptide (I = L), all GENCODE CDSs and lincorfz ORFs whose translation contains it are
found; the encoding codons are extracted and every distinct nucleotide sequence becomes one read
(name: <peptide>|<variant>|<n_sources>). A full 3' adapter is appended so that the adapter
trimming step of the alignment workflow removes exactly the adapter and keeps the sequence.
Peptides are assigned to read libraries by the PSMs (Spectrum File) and, for DIA data, by the DIA-NN
precursor matrix (1% FDR); one FASTQ is written per library. A library is given as a substring of the
spectrum file / run column names, or as <name>=<regex> to pool several runs.

Outputs (in --outdir):
  <run>_R1.fastq.gz                  reads for the peptides detected in that run
  peptide_sources.tsv                peptide, source (transcript / ORF), gene, offset, nucleotide sequence
  peptides.tsv                       peptide, runs, number of sources and distinct encodings
  samplesheet.tsv                    sample, fq1, expected_size (input of the alignment workflow)
"""
import argparse
import gzip
import os
import re
from collections import defaultdict

import pandas as pd
import pysam

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
COMP = str.maketrans("ACGTN", "TGCAN")
IL = str.maketrans("I", "L")
K = 7  # minimal peptide length = prefix length of the index


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def revcomp(s):
    return s.translate(COMP)[::-1]


def gencode_cds(gtf, fa):
    """Yield (transcript_id, gene_name, nt) of every CDS in the GTF (phase-corrected)."""
    cds = defaultdict(list)
    info = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "CDS":
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            gname = re.search(r'gene_name "([^"]+)"', f[8])
            cds[tid].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] != "." else 0))
            info[tid] = (f[0], f[6], gname.group(1) if gname else tid)
    for tid, segs in cds.items():
        chrom, strand, gname = info[tid]
        if chrom not in fa.references:
            continue
        segs.sort(reverse=(strand == "-"))
        parts = [fa.fetch(chrom, s, e).upper() for s, e, _ in segs]
        nt = "".join(revcomp(p) for p in parts) if strand == "-" else "".join(parts)
        yield tid, gname, nt[segs[0][2]:]


def lincorf_nt(orfs_tsv, fa):
    o = pd.read_csv(orfs_tsv, sep="\t", usecols=["orf_label", "locus_name", "chrom", "strand", "orf_blocks"],
                    low_memory=False)
    for r in o.itertuples():
        blocks = [tuple(map(int, b.split("-"))) for b in re.split(r"[,;]", str(r.orf_blocks)) if b]
        parts = [fa.fetch(r.chrom, s, e).upper() for s, e in sorted(blocks)]
        nt = "".join(parts)
        yield r.orf_label, r.locus_name, (revcomp(nt) if r.strand == "-" else nt)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psm", required=True, nargs="+", help="FragPipe psm.tsv file(s)")
    ap.add_argument("--pr-matrix", help="DIA-NN report.pr_matrix.tsv (DIA only)")
    ap.add_argument("--runs", required=True,
                    help="comma-separated libraries: run name substring, or <name>=<regex> over spectrum files/runs")
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (CDS features)")
    ap.add_argument("--genome", required=True, help="genome FASTA (indexed), same assembly as the GTF")
    ap.add_argument("--orfs", required=True, help="lincorfz lincORFs.ranked.tsv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    libs = {}
    for item in args.runs.split(","):
        name, _, rx = item.partition("=")
        libs[name] = re.compile(rx if rx else re.escape(name))
    runs = list(libs)

    # peptides per library
    pep_runs = defaultdict(set)
    if args.pr_matrix:
        pr = pd.read_csv(args.pr_matrix, sep="\t", low_memory=False)
        for run, rx in libs.items():
            cols = [c for c in pr.columns if rx.search(c)]
            for p in pr.loc[pr[cols].notna().any(axis=1), "Stripped.Sequence"]:
                pep_runs[p].add(run)
    psm = pd.concat([pd.read_csv(f, sep="\t", usecols=["Spectrum File", "Peptide"], low_memory=False)
                     for f in args.psm], ignore_index=True)
    for run, rx in libs.items():
        for p in psm.loc[psm["Spectrum File"].astype(str).str.contains(rx), "Peptide"]:
            pep_runs[p].add(run)
    peps = sorted(p for p in pep_runs if len(p) >= K)
    by_prefix = defaultdict(list)
    for p in peps:
        by_prefix[p.translate(IL)[:K]].append(p)
    print(len(peps), "peptides")

    # scan all coding sequences once with a K-mer window, verify full peptides on prefix hits
    fa = pysam.FastaFile(args.genome)
    hits = defaultdict(list)    # peptide -> [(source, gene, type, aa_offset, nt)]
    for src_type, gen in [("GENCODE_CDS", gencode_cds(args.gtf, fa)), ("lincORF", lincorf_nt(args.orfs, fa))]:
        n = 0
        for sid, gname, nt in gen:
            n += 1
            aa = translate(nt).translate(IL)
            for i in range(len(aa) - K + 1):
                cand = by_prefix.get(aa[i:i + K])
                if not cand:
                    continue
                for p in cand:
                    q = p.translate(IL)
                    if aa.startswith(q, i):
                        hits[p].append((sid, gname, src_type, i, nt[3 * i:3 * (i + len(p))]))
        print(src_type, n, "sequences scanned")

    rows, pep_rows = [], []
    reads = defaultdict(list)   # run -> fastq records
    for p in peps:
        h = hits.get(p, [])
        variants = {}
        for sid, gname, st, off, nt in h:
            variants.setdefault(nt, []).append((sid, gname, st))
            rows.append((p, sid, gname, st, off, nt))
        pep_rows.append((p, ",".join(sorted(pep_runs[p])), len(h), len({x[1] for x in h}), len(variants),
                         ",".join(sorted({x[2] for x in h}))))
        for k, (nt, srcs) in enumerate(sorted(variants.items())):
            name = f"{p}|v{k + 1}|{len(srcs)}"
            seq = nt + ADAPTER
            for run in pep_runs[p]:
                reads[run].append(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
    pd.DataFrame(rows, columns=["peptide", "source_id", "gene_name", "source_type", "aa_offset", "nt_seq"]) \
        .to_csv(f"{args.outdir}/peptide_sources.tsv", sep="\t", index=False)
    pd.DataFrame(pep_rows, columns=["peptide", "runs", "n_sources", "n_genes", "n_encodings", "source_types"]) \
        .to_csv(f"{args.outdir}/peptides.tsv", sep="\t", index=False)
    sheet = []
    for run in runs:
        path = os.path.abspath(f"{args.outdir}/{run}_R1.fastq.gz")
        with gzip.open(path, "wt") as fh:
            fh.writelines(reads[run])
        sheet.append((run, path, os.path.getsize(path)))
        print(run, len(reads[run]), "reads")
    with open(f"{args.outdir}/samplesheet.tsv", "w") as fh:
        fh.write("sample\tfq1\texpected_size\n")
        fh.writelines(f"{s}\t{p}\t{z}\n" for s, p, z in sheet)
    n_src = sum(1 for r in pep_rows if r[2] > 0)
    print(f"{n_src}/{len(pep_rows)} peptides with at least one encoding sequence")


if __name__ == "__main__":
    main()
