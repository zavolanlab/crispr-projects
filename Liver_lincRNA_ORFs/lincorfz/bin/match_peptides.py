#!/usr/bin/env python
"""
Optional: match identified peptides (e.g. from non-tryptic small-peptidome searches) to ORFs.

Peptide tables can be any TSV/CSV with a column named one of
'Peptide', 'peptide', 'Sequence', 'sequence', 'Stripped.Sequence', 'Peptide Sequence';
modifications in brackets/parentheses and non-letters are stripped. Matching treats I and L
as equivalent. Peptides that also occur in the canonical proteome (FASTA, e.g. translated
GENCODE CDSs) are discarded as non-unique. Peptides shorter than --min-length are ignored.

Output: peptides.orfs.tsv (orf_id, n_unique_peptides, peptides, datasets)
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

PEP_COLS = ["Peptide", "peptide", "Sequence", "sequence", "Stripped.Sequence", "Peptide Sequence",
            "Stripped Sequence"]


def norm(p):
    p = re.sub(r"\[.*?\]|\(.*?\)", "", str(p))
    return re.sub(r"[^A-Z]", "", p.upper()).replace("I", "L")


def read_fasta(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orfs", required=True)
    ap.add_argument("--peptide-tables", nargs="+", required=True)
    ap.add_argument("--proteome-fasta", required=True)
    ap.add_argument("--min-length", type=int, default=7)
    args = ap.parse_args()

    peps = {}
    for path in args.peptide_tables:
        sep = "," if path.endswith(".csv") else "\t"
        df = pd.read_csv(path, sep=sep, low_memory=False)
        col = next((c for c in PEP_COLS if c in df.columns), None)
        if col is None:
            raise SystemExit(f"{path}: no peptide column among {PEP_COLS}")
        for p in df[col].dropna().map(norm):
            if len(p) >= args.min_length:
                peps.setdefault(p, set()).add(os.path.basename(path))
    proteome = "|".join(s.replace("I", "L") for s in read_fasta(args.proteome_fasta))
    unique = {p: d for p, d in peps.items() if p not in proteome}
    print(f"{len(peps)} peptides, {len(unique)} absent from the canonical proteome")

    orfs = pd.read_csv(args.orfs, sep="\t", usecols=["orf_id", "aa_seq"])
    # one concatenated string; occurrences are mapped back to ORFs through the offsets
    seqs = [str(s).replace("I", "L") for s in orfs["aa_seq"]]
    starts = np.cumsum([0] + [len(s) + 1 for s in seqs[:-1]])
    concat = "|".join(seqs)
    hits = {}
    for p in unique:
        i = concat.find(p)
        while i >= 0:
            hits.setdefault(int(np.searchsorted(starts, i, side="right") - 1), set()).add(p)
            i = concat.find(p, i + 1)
    rows = []
    for idx, ps in sorted(hits.items()):
        rows.append(dict(orf_id=orfs["orf_id"].iat[idx], n_unique_peptides=len(ps), peptides=",".join(sorted(ps)),
                         datasets=",".join(sorted(set().union(*[unique[p] for p in ps])))))
    pd.DataFrame(rows, columns=["orf_id", "n_unique_peptides", "peptides", "datasets"]).to_csv(
        "peptides.orfs.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
