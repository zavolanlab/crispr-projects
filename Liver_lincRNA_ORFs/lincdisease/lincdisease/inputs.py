"""
Inputs of lincdisease (LOGIC.md, "Scope and inputs"): the candidate ORF table, taken as is, and the
background of all ORFs / loci of the same lincorfz run (results/final/lincORFs.ranked.tsv).
Candidates are part of the background (percentiles are computed against all expressed lincRNA ORFs).
"""
import os

import pandas as pd


def read_table(path):
    """lincorfz tables are tab-separated whatever the extension (the curated list ends in .tsv.csv)."""
    return pd.read_csv(path, sep="\t", low_memory=False)


def load(results_dir, candidates):
    """
    Background ORFs with an `is_candidate` column. `candidates` is a path (any lincorfz-format table)
    or a list of orf_ids. ORFs of the candidate table missing from the background are an error.
    """
    bg = read_table(os.path.join(results_dir, "final", "lincORFs.ranked.tsv"))
    if isinstance(candidates, str):
        cand_ids = read_table(candidates)["orf_id"]
    else:
        cand_ids = pd.Series(list(candidates))
    missing = set(cand_ids) - set(bg["orf_id"])
    if missing:
        raise ValueError(f"{len(missing)} candidate ORFs not in the background, e.g. {sorted(missing)[:3]}")
    bg = bg.copy()
    bg["is_candidate"] = bg["orf_id"].isin(set(cand_ids))
    return bg


def loci(orfs):
    """One row per locus: the locus columns of its best-ranked ORF and whether it carries a candidate."""
    first = orfs.sort_values("rank").drop_duplicates("locus_id")
    cols = ["locus_id", "locus_name", "gene_name_source", "name_source_gene_ids", "chrom", "strand",
            "lncRNA_gene_names", "lncRNA_gene_ids", "small_RNA_names", "antisense_coding_exon_overlap",
            "log2TPM_HepG2", "log2TPM_Huh7", "TPM_HepG2", "TPM_Huh7"]
    out = first[cols].set_index("locus_id")
    out["is_candidate"] = orfs.groupby("locus_id")["is_candidate"].any()
    return out
