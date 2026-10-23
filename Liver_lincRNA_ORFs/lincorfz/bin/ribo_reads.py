"""
Read-level helpers for genomic Ribo-seq alignments (pysam AlignedSegment).

Reads are stratified by (aligned length, 5' soft-clip class). The 5' soft clip
captures untemplated nucleotides added during library preparation, which shift the
aligned 5' end relative to the true footprint 5' end; P-site offsets are therefore
calibrated per class and always counted from the *aligned* 5' end.
"""
import gzip

import pandas as pd

BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CEQUAL, BAM_CDIFF = 0, 1, 2, 3, 4, 7, 8
REF_CONSUMING = {BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP, BAM_CEQUAL, BAM_CDIFF}


def clip5_of(read, max_clip5):
    cig = read.cigartuples
    op, ln = cig[-1] if read.is_reverse else cig[0]
    c = ln if op == BAM_CSOFT_CLIP else 0
    return min(c, max_clip5)


def introns_of_read(read):
    """0-based half-open introns (N operations) of an alignment."""
    out = []
    pos = read.reference_start
    for op, ln in read.cigartuples:
        if op == BAM_CREF_SKIP:
            out.append((pos, pos + ln))
        if op in REF_CONSUMING:
            pos += ln
    return out


def has_skip(read):
    return any(op == BAM_CREF_SKIP for op, _ in read.cigartuples)


def five_prime_genomic(read):
    """0-based genomic coordinate of the aligned 5' end (strand-aware)."""
    return read.reference_end - 1 if read.is_reverse else read.reference_start


def psite_genomic(read, offset):
    """
    Genomic 0-based position located `offset` aligned reference nucleotides downstream
    of the aligned 5' end, walking along the alignment blocks (splice-aware).
    Returns None if the offset falls beyond the aligned part of the read.
    """
    blocks = read.get_blocks()
    if read.is_reverse:
        rem = offset
        for s, e in reversed(blocks):
            ln = e - s
            if rem < ln:
                return e - 1 - rem
            rem -= ln
        return None
    rem = offset
    for s, e in blocks:
        ln = e - s
        if rem < ln:
            return s + rem
        rem -= ln
    return None


def load_whitelist(path, chroms=None):
    wl = pd.read_csv(path, sep="\t")
    if chroms is not None:
        wl = wl.loc[wl["chrom"].isin(chroms)]
    out = {}
    for c, g in wl.groupby("chrom"):
        out[c] = set(zip(g["start0"].astype(int), g["end1"].astype(int)))
    return out


def junctions_ok(read, chrom_whitelist):
    return all(j in chrom_whitelist for j in introns_of_read(read))


def load_offsets(path):
    """dict (aligned_len, clip5) -> offset for classes that passed calibration."""
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["pass"].astype(bool)]
    return {(int(r.aligned_len), int(r.clip5)): int(r.offset) for r in df.itertuples()}
