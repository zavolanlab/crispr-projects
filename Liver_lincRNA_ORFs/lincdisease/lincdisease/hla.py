"""
Predicted HLA class I presentation of lincORF peptides (actionability score, LOGIC.md Step 3) with MHCflurry 2
(Class1PresentationPredictor: binding affinity + antigen processing, O'Donnell et al. 2020).

All 8-11-mers of each ORF peptide are scored for three sets of alleles:
* HepG2 (Cellosaurus CVCL_0027, PubMed 25960936): A*02:01, A*24:02, B*35:14, B*51:01, C*04:01, C*16:02;
* Huh7 (Cellosaurus CVCL_0336, CLS typing; homozygous): A*11:01, B*54:01, C*01:02;
* common alleles: representatives of the main HLA-I supertypes plus alleles frequent in East Asian HCC
  populations, each scored as its own one-allele genotype (COMMON below).
A peptide is a predicted binder if its presentation percentile is <= 2 (MHCflurry's usual cut-off; percentile ranks
are calibrated per allele). Peptides that also occur in an annotated protein (GENCODE v42 protein-coding translations,
UniProt human reference proteome with isoforms) are self-peptides, to which T cells are tolerant, and are left out.
Per ORF and allele set: number of distinct non-self binder peptides and the best percentile of the non-self peptides;
hla_n_self_binders counts the binders left out.

MHCflurry runs in its own environment; `predict` is executed there via `python -m lincdisease.hla` (see main()).
"""
import gzip
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

GENOTYPES = {
    "HepG2": ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*35:14", "HLA-B*51:01", "HLA-C*04:01", "HLA-C*16:02"],
    "Huh7": ["HLA-A*11:01", "HLA-B*54:01", "HLA-C*01:02"],
}
COMMON = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02", "HLA-A*33:03",
          "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*27:05", "HLA-B*40:01", "HLA-B*44:02",
          "HLA-B*46:01", "HLA-B*58:01", "HLA-C*07:01", "HLA-C*07:02"]
BINDER_PERCENTILE = 2.0
LENGTHS = (8, 9, 10, 11)
# actionability control: GPC3 144-152, HLA-A*02:01 peptide vaccine epitope
CONTROL = ("GPC3_144-152", "FVGEFFTDV", "HLA-A*02:01")


def read_proteins(paths):
    """Protein sequences of (gzipped) FASTA files, deduplicated."""
    seqs, cur = set(), []
    for path in paths:
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        seqs.add("".join(cur))
                    cur = []
                else:
                    cur.append(line.strip())
        if cur:
            seqs.add("".join(cur))
            cur = []
    return seqs


def self_peptides(seqs, proteins, lengths=LENGTHS):
    """The 8-11-mers of `seqs` that occur in any of `proteins` (exact match), found by one scan of the proteins over
    the k-mer prefixes of the query peptides."""
    k = min(lengths)
    by_prefix = defaultdict(set)
    for s in seqs:
        for n in lengths:
            for i in range(len(s) - n + 1):
                by_prefix[s[i:i + k]].add(s[i:i + n])
    found = set()
    for prot in proteins:
        for i in range(len(prot) - k + 1):
            cand = by_prefix.get(prot[i:i + k])
            if cand:
                found.update(p for p in cand if prot.startswith(p, i))
    return found


def summarize(pred, genotype_of_sample, threshold=BINDER_PERCENTILE, self_peps=frozenset()):
    """
    MHCflurry `predict_sequences(result='all')` rows (sequence_name, sample_name, peptide, presentation_percentile)
    -> per sequence: <set>_n_binders and <set>_best_pct for every allele set, plus any_n_binders (distinct binder
    peptides over all sets) and n_self_binders (binders in `self_peps`, left out of the other columns).
    """
    p = pred.assign(set=pred["sample_name"].map(genotype_of_sample))
    is_self = p["peptide"].isin(self_peps)
    n_self = p[is_self & (p["presentation_percentile"] <= threshold)].groupby("sequence_name")["peptide"].nunique()
    p = p[~is_self]
    best = p.groupby(["sequence_name", "set"])["presentation_percentile"].min().unstack()
    binders = p[p["presentation_percentile"] <= threshold]
    n = binders.groupby(["sequence_name", "set"])["peptide"].nunique().unstack()
    n = n.reindex(index=best.index, columns=best.columns).fillna(0).astype(int)
    out = pd.concat([n.add_suffix("_n_binders"), best.add_suffix("_best_pct")], axis=1)
    out["any_n_binders"] = binders.groupby("sequence_name")["peptide"].nunique().reindex(out.index).fillna(0).astype(int)
    out["n_self_binders"] = n_self.reindex(out.index).fillna(0).astype(int)
    out.columns = ["hla_" + c for c in out.columns]
    return out


def samples():
    """MHCflurry sample -> alleles, and sample -> allele set name."""
    alleles = dict(GENOTYPES)
    alleles.update({f"common:{a}": [a] for a in COMMON})
    return alleles, {s: ("common" if s.startswith("common:") else s) for s in alleles}


def predict(seqs, chunk=500, control_path=None, protein_fastas=()):
    """Run MHCflurry on Series key -> aa sequence (only in the MHCflurry environment); the GPC3 control prediction
    is written to `control_path` if given. Peptides found in `protein_fastas` are self-peptides (summarize)."""
    from mhcflurry import Class1PresentationPredictor

    predictor = Class1PresentationPredictor.load()
    alleles, set_of = samples()
    seqs = seqs[seqs.str.len() >= min(LENGTHS)]
    self_peps = self_peptides(seqs, read_proteins(protein_fastas)) if protein_fastas else frozenset()
    print(f"self-peptides among the ORF 8-11-mers: {len(self_peps)}", flush=True)
    parts = []
    keys = list(seqs.index)
    for i in range(0, len(keys), chunk):
        sub = seqs[keys[i:i + chunk]].to_dict()
        r = predictor.predict_sequences(sub, alleles=alleles, result="all", peptide_lengths=LENGTHS, verbose=0)
        parts.append(summarize(r, set_of, self_peps=self_peps))
        print(f"{min(i + chunk, len(keys))}/{len(keys)}", flush=True)
    ctrl = predictor.predict([CONTROL[1]], alleles=[CONTROL[2]], verbose=0)
    print(f"control {CONTROL[0]} {CONTROL[1]} {CONTROL[2]}: affinity {ctrl['affinity'].iloc[0]:.1f} nM, "
          f"presentation percentile {ctrl['presentation_percentile'].iloc[0]:.3f}", flush=True)
    if control_path:
        ctrl[["peptide", "affinity", "presentation_percentile"]].assign(control=CONTROL[0], allele=CONTROL[2]).to_csv(
            control_path, sep="\t", index=False)
    return pd.concat(parts)


def run(orfs, python, out_path, env=None, protein_fastas=()):
    """Write the ORF peptides, run `python -m lincdisease.hla` in the MHCflurry environment, read the result."""
    out_path = os.path.abspath(out_path)     # the subprocess runs in the package directory
    seq_path = out_path + ".input.tsv"
    seqs = orfs.set_index("orf_id")["aa_seq"].fillna("").str.rstrip("*")
    seqs.rename("aa_seq").to_csv(seq_path, sep="\t")
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    subprocess.run([python, "-m", "lincdisease.hla", seq_path, out_path, *map(os.path.abspath, protein_fastas)],
                   check=True, cwd=pkg_root,
                   env={**os.environ, **(env or {})})
    return pd.read_csv(out_path, sep="\t", index_col=0)


def finalize(res, orf_ids):
    """One row per ORF; ORFs shorter than 8 aa have no peptides (NaN percentiles, 0 binders)."""
    res = res.reindex(orf_ids)
    n_cols = [c for c in res.columns if c.endswith("_n_binders")]
    res[n_cols] = res[n_cols].fillna(0).astype(int)
    return res


def main():
    seqs = pd.read_csv(sys.argv[1], sep="\t", index_col=0)["aa_seq"].fillna("")
    control_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[2])), "hla.control.tsv")
    res = finalize(predict(seqs, control_path=control_path, protein_fastas=sys.argv[3:]), seqs.index)
    res.to_csv(sys.argv[2], sep="\t", index_label="orf_id")


if __name__ == "__main__":
    main()
