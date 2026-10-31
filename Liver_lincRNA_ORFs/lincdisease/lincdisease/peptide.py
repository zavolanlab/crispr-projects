"""
M1 (encoded peptide) sequence and annotation features of every ORF (LOGIC.md, Step 4).

* Independent translation evidence: the ORF has the same stop codon (chrom, strand, genomic stop end) as an ORF
  of the GENCODE Ribo-seq ORF catalog (phase 2 comprehensive set, GENCODE v45; phase 1 ORFs with peptide support
  reported separately). Both catalogs and lincorfz include the stop codon in the ORF coordinates.
* Conservation contrast: mean PhyloP100way over the ORF minus the mean over the exonic positions of its transcript
  outside the ORF (a translated and functional peptide should be more constrained than the rest of the RNA).
* Transmembrane segments: Kyte-Doolittle hydropathy, 19-residue windows with mean >= 1.6 (Kyte & Doolittle 1982);
  overlapping windows form one segment. Open, transparent, and adequate for the short ORFs here.
* Signal-peptide-like N-terminus (heuristic, as open predictors such as DeepSig no longer install and SignalP is
  licensed): a hydrophobic h-region (8-residue window, mean hydropathy >= 2.0, no D/E/K/R) starting within
  residues 2-20, followed 3-12 residues downstream by a small residue pair at -3 and -1 (A/G/S/C/T; von Heijne's
  -3,-1 rule) at a cleavage position <= 40.
* Pfam domains: `hmmsearch --cut_ga` against Pfam-A (domain-level gathering thresholds).
* Disorder: metapredict (per-residue score >= 0.5 is disordered); fraction of disordered residues.
* Length in amino acids.

External tools run in their own conda environment (`hmmsearch`, a Python with metapredict), called by path.
"""
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
TM_WINDOW, TM_THRESHOLD = 19, 1.6
SMALL = set("AGSCT")
PHYLOP_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"


# ---------------------------------------------------------------- coordinates

def parse_blocks(s):
    """lincorfz `orf_blocks` ('start0-end1;start0-end1') -> list of (start0, end1)."""
    return [tuple(int(x) for x in b.split("-")) for b in str(s).split(";") if b]


def read_bed12(path):
    """BED12 -> chrom/strand/name/blocks (list of (start0, end1))."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12 or ln.startswith(("#", "track")):
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].strip(",").split(",")]
            starts = [int(x) for x in f[11].strip(",").split(",")]
            rows.append((f[0], f[5], f[3], [(start + s, start + s + n) for s, n in zip(starts, sizes)]))
    return pd.DataFrame(rows, columns=["chrom", "strand", "name", "blocks"])


def stop_end(blocks, strand):
    """Genomic end of the stop codon: 3'-most coordinate of the ORF on its strand."""
    return blocks[-1][1] if strand == "+" else blocks[0][0]


def same_stop(orfs, catalog):
    """Per ORF: ';'-joined names of catalog ORFs ending at the same stop codon ('' if none)."""
    key = {}
    for c in catalog.itertuples(index=False):
        key.setdefault((c.chrom, c.strand, stop_end(c.blocks, c.strand)), []).append(c.name.split(":")[0])
    out = [";".join(key.get((o.chrom, o.strand, stop_end(parse_blocks(o.orf_blocks), o.strand)), []))
           for o in orfs.itertuples(index=False)]
    return pd.Series(out, index=orfs.index)


# ---------------------------------------------------------------- sequence features

def hydropathy(seq, window):
    """Mean Kyte-Doolittle hydropathy of every window (length len(seq) - window + 1)."""
    v = np.array([KD.get(a, 0.0) for a in seq])
    if len(v) < window:
        return np.array([])
    c = np.concatenate([[0], np.cumsum(v)])
    return (c[window:] - c[:-window]) / window


def tm_segments(seq, window=TM_WINDOW, threshold=TM_THRESHOLD):
    """Number of predicted TM segments: runs of overlapping windows with mean hydropathy >= threshold."""
    h = hydropathy(seq, window) >= threshold
    if not h.any():
        return 0
    return int(np.sum(h[1:] & ~h[:-1]) + h[0])


def signal_like(seq, h_window=8, h_threshold=2.0, h_start_max=20, max_cleavage=40):
    """Heuristic signal-peptide-like N-terminus (see module docstring)."""
    s = seq.rstrip("*")
    h = hydropathy(s[:h_start_max + h_window], h_window)
    for i in np.flatnonzero(h >= h_threshold):
        if i < 1 or set(s[i:i + h_window]) & set("DEKR"):
            continue
        h_end = i + h_window                 # 0-based index after the h-region
        for c in range(h_end + 3, min(h_end + 13, max_cleavage + 1, len(s))):
            # cleavage between residues c-1 and c (0-based): -1 = s[c-1], -3 = s[c-3]
            if s[c - 1] in SMALL and s[c - 3] in SMALL:
                return True
    return False


def sequence_features(orfs):
    """Length, TM segments, max 19-window hydropathy and signal-peptide-like N-terminus from `aa_seq`."""
    seq = orfs["aa_seq"].fillna("").str.rstrip("*")
    return pd.DataFrame({
        "pep_length": seq.str.len(),
        "pep_tm_segments": seq.map(tm_segments),
        "pep_max_hydropathy19": seq.map(lambda s: float(hydropathy(s, TM_WINDOW).max()) if len(s) >= TM_WINDOW
                                        else np.nan),
        "pep_signal_like": seq.map(signal_like),
    }, index=orfs.index)


# ---------------------------------------------------------------- external tools

def write_fasta(seqs, path):
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")


DISORDER_SCRIPT = """
import sys, metapredict as meta
ids, seqs = [], []
for ln in open(sys.argv[1]):
    if ln.startswith(">"): ids.append(ln[1:].strip())
    else: seqs.append(ln.strip())
with open(sys.argv[2], "w") as out:
    out.write("key\\tmean_disorder\\tfrac_disordered\\n")
    for k, s in zip(ids, seqs):
        d = meta.predict_disorder(s)
        out.write(f"{k}\\t{sum(d) / len(d):.4f}\\t{sum(x >= 0.5 for x in d) / len(d):.4f}\\n")
"""


def disorder(seqs, python):
    """metapredict per sequence (Series key -> aa) -> DataFrame mean_disorder, frac_disordered (index = key)."""
    with tempfile.TemporaryDirectory() as d:
        fa, out = os.path.join(d, "in.fa"), os.path.join(d, "out.tsv")
        write_fasta(seqs, fa)
        subprocess.run([python, "-c", DISORDER_SCRIPT, fa, out], check=True,
                       env={**os.environ, "MPLCONFIGDIR": d})
        return pd.read_csv(out, sep="\t", index_col=0)


def read_domtblout(path):
    """hmmsearch --domtblout -> target (sequence), pfam name, accession, i-Evalue, envelope."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split()
            rows.append((f[0], f[3], f[4].split(".")[0], float(f[12]), int(f[19]), int(f[20])))
    return pd.DataFrame(rows, columns=["key", "pfam_name", "pfam_acc", "i_evalue", "env_from", "env_to"])


def pfam(seqs, hmmsearch, hmm, cpu=4):
    """Pfam-A domains (gathering thresholds) per sequence key: number of domains and ';'-joined names."""
    with tempfile.TemporaryDirectory() as d:
        fa, dom = os.path.join(d, "in.fa"), os.path.join(d, "dom.tsv")
        write_fasta(seqs, fa)
        subprocess.run([hmmsearch, "--cut_ga", "--noali", "--cpu", str(cpu), "--domtblout", dom, "-o", os.devnull,
                        hmm, fa], check=True)
        hits = read_domtblout(dom)
    g = hits.groupby("key")
    return pd.DataFrame({"pfam_n_domains": g.size(),
                         "pfam_domains": g["pfam_name"].agg(lambda s: ";".join(sorted(set(s))))})


# ---------------------------------------------------------------- conservation contrast

def _positions(blocks):
    return np.concatenate([np.arange(s, e) for s, e in blocks]) if blocks else np.array([], dtype=int)


def conservation_contrast(orfs, tx_exons, bigwig=PHYLOP_URL, threads=8):
    """
    Per ORF: mean PhyloP over the ORF, over the exonic positions of its transcript outside the ORF, and their
    difference. One bigWig query per locus span (works on a remote URL).
    """
    import pyBigWig

    ex = tx_exons[tx_exons["transcript_id"].isin(set(orfs["transcript_id"]))]
    tx_blocks = {t: list(zip(g["start0"], g["end1"])) for t, g in ex.groupby("transcript_id")}
    groups = list(orfs.groupby("locus_id"))

    def one(item):
        _, g = item
        bw = pyBigWig.open(bigwig)
        try:
            chrom = g["chrom"].iloc[0]
            spans = [b for t in g["transcript_id"] for b in tx_blocks.get(t, [])]
            spans += [b for s in g["orf_blocks"] for b in parse_blocks(s)]
            lo, hi = min(s for s, _ in spans), max(e for _, e in spans)
            vals = np.array(bw.values(chrom, lo, hi), dtype=float)
        finally:
            bw.close()
        res = []
        for o in g.itertuples():
            orf_pos = _positions(parse_blocks(o.orf_blocks))
            rest = np.setdiff1d(_positions(tx_blocks.get(o.transcript_id, [])), orf_pos)
            a, b = vals[orf_pos - lo], vals[rest - lo]
            ma = np.nanmean(a) if np.isfinite(a).any() else np.nan
            mb = np.nanmean(b) if np.isfinite(b).sum() >= 30 else np.nan
            res.append((o.Index, ma, mb, len(rest)))
        return res

    with ThreadPoolExecutor(threads) as pool:
        rows = [r for part in pool.map(one, groups) for r in part]
    out = pd.DataFrame(rows, columns=["idx", "phylop_orf", "phylop_nonorf_exons", "nonorf_exon_nt"]).set_index("idx")
    out["phylop_orf_minus_nonorf"] = out["phylop_orf"] - out["phylop_nonorf_exons"]
    return out.reindex(orfs.index)
