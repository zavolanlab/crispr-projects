"""
Shared helpers for the lincorfz pipeline scripts.

Coordinate conventions used throughout the pipeline:
  * exons are stored as 0-based half-open genomic intervals (start0, end1),
    sorted by genomic start, independently of strand;
  * "transcript coordinates" are 0-based offsets from the transcript 5' end
    (i.e. they run along the exon chain in the direction of transcription).
"""
import gzip
import re

import numpy as np
import pandas as pd

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ----------------------------------------------------------------------------------
# generic
# ----------------------------------------------------------------------------------
def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_chromosomes(arg, chrom_sizes):
    """'all' -> main chromosomes (chr1-22, X, Y, M) present in chrom_sizes; else comma list."""
    sizes = pd.read_csv(chrom_sizes, sep="\t", header=None, names=["chrom", "size"])
    if arg in (None, "", "all"):
        main = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y", "M"]]
        return [c for c in main if c in set(sizes["chrom"])]
    return [c.strip() for c in str(arg).replace(" ", ",").split(",") if c.strip()]


def bigwig_chrom_name(bw, chrom):
    """Return the chromosome name used inside a bigwig for a 'chr'-style name (or None)."""
    chroms = bw.chroms()
    candidates = [chrom]
    if chrom.startswith("chr"):
        bare = chrom[3:]
        candidates += ["MT" if bare == "M" else bare]
    else:
        candidates += ["chrM" if chrom == "MT" else "chr" + chrom]
    for c in candidates:
        if c in chroms:
            return c
    return None


# ----------------------------------------------------------------------------------
# GTF
# ----------------------------------------------------------------------------------
def parse_attributes(attr_string, keys=None):
    d = {}
    for k, v in ATTR_RE.findall(attr_string):
        if keys is None or k in keys:
            # GTF allows repeated keys (e.g. tag); keep them comma-joined
            d[k] = d[k] + "," + v if k in d else v
    return d


def read_gtf(path, features=None, keys=("gene_id", "transcript_id", "gene_name", "gene_type",
                                        "transcript_type", "tag", "ref_gene_id", "ref_gene_name",
                                        "cmp_ref", "class_code", "TPM", "FPKM", "cov")):
    """Read a GTF into a DataFrame with 0-based 'start0' and 1-based 'end1' columns."""
    rows = []
    keys = set(keys) if keys is not None else None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if features is not None and f[2] not in features:
                continue
            rec = {"chrom": f[0], "source": f[1], "feature": f[2], "start0": int(f[3]) - 1,
                   "end1": int(f[4]), "score": f[5], "strand": f[6], "frame": f[7]}
            rec.update(parse_attributes(f[8], keys))
            rows.append(rec)
    df = pd.DataFrame(rows)
    return df


def exon_chains(exons_df, id_col="transcript_id"):
    """dict id -> (chrom, strand, [(start0, end1), ...] sorted by genomic start)."""
    chains = {}
    for tid, g in exons_df.groupby(id_col, sort=False):
        g = g.sort_values("start0")
        chains[tid] = (g["chrom"].iat[0], g["strand"].iat[0],
                       list(zip(g["start0"].astype(int), g["end1"].astype(int))))
    return chains


def introns_of_chain(exons):
    """0-based half-open introns (start0, end1) between consecutive exons."""
    return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)
            if exons[i + 1][0] > exons[i][1]]


# ----------------------------------------------------------------------------------
# transcript <-> genome mapping
# ----------------------------------------------------------------------------------
def tx_genomic_positions(exons, strand):
    """numpy array: genomic 0-based position of each transcript nucleotide (5'->3')."""
    pos = np.concatenate([np.arange(s, e) for s, e in exons]) if exons else np.array([], dtype=int)
    return pos[::-1] if strand == "-" else pos


class TxMapper:
    """Fast genomic position -> transcript coordinate lookup for one transcript."""

    def __init__(self, exons, strand):
        self.exons = exons
        self.strand = strand
        self.starts = np.array([s for s, _ in exons])
        self.ends = np.array([e for _, e in exons])
        lens = self.ends - self.starts
        self.length = int(lens.sum())
        # cumulative length before each exon, in genomic order
        self.cum = np.concatenate([[0], np.cumsum(lens)[:-1]])

    def to_tx(self, gpos):
        """genomic 0-based position -> transcript coordinate, or -1 if not exonic."""
        i = np.searchsorted(self.starts, gpos, side="right") - 1
        if i < 0 or gpos >= self.ends[i]:
            return -1
        g_off = self.cum[i] + (gpos - self.starts[i])
        return int(self.length - 1 - g_off) if self.strand == "-" else int(g_off)


def blocks_from_positions(gpos):
    """Collapse sorted genomic positions into (start0, end1) blocks."""
    gpos = np.sort(np.asarray(gpos))
    if len(gpos) == 0:
        return []
    breaks = np.flatnonzero(np.diff(gpos) != 1)
    starts = np.concatenate([[gpos[0]], gpos[breaks + 1]])
    ends = np.concatenate([gpos[breaks] + 1, [gpos[-1] + 1]])
    return list(zip(starts.astype(int), ends.astype(int)))


def spliced_sequence(fasta, chrom, exons, strand):
    seq = "".join(fasta.fetch(chrom, s, e) for s, e in exons).upper()
    return revcomp(seq) if strand == "-" else seq


def bigwig_along_tx(bw, chrom_in_bw, exons, strand):
    """Values of a bigwig along the exon chain, in transcript orientation (NaN -> 0)."""
    if chrom_in_bw is None:
        return np.zeros(sum(e - s for s, e in exons))
    vals = [np.nan_to_num(bw.values(chrom_in_bw, s, e, numpy=True), nan=0.0) for s, e in exons]
    v = np.concatenate(vals) if vals else np.array([])
    return v[::-1] if strand == "-" else v


# ----------------------------------------------------------------------------------
# genetic code
# ----------------------------------------------------------------------------------
_bases = "TCAG"
_codons = [a + b + c for a in _bases for b in _bases for c in _bases]
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(_codons, _aas))
STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - len(seq) % 3, 3))
