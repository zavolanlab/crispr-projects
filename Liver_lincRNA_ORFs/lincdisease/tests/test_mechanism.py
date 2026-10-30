import numpy as np
import pandas as pd
import pytest

from lincdisease import mechanism


def _tx(strand):
    # two exons: 0-100 and 200-300
    return pd.DataFrame({"transcript_id": ["T", "T"], "chrom": "c1", "start0": [0, 200], "end1": [100, 300],
                         "strand": strand})


def test_stop_to_last_junction_plus_and_minus():
    orfs = pd.DataFrame({"transcript_id": ["T", "T"], "orf_start0": [10, 220], "orf_end1": [40, 260],
                         "strand": ["+", "+"]})
    d = mechanism.stop_to_last_junction(orfs, _tx("+"))
    assert d.tolist() == [60, 0]          # 40..100 is exonic upstream of the junction at 200
    orfs_m = pd.DataFrame({"transcript_id": ["T"], "orf_start0": [250], "orf_end1": [290], "strand": ["-"]})
    assert mechanism.stop_to_last_junction(orfs_m, _tx("-")).tolist() == [50]   # 200..250


def test_neighbour_position_divergent_and_antisense():
    loci = pd.DataFrame({"gene_id": ["div", "anti", "far"], "chrom": "c1", "start0": [10_000, 50_000, 90_000],
                         "end1": [12_000, 52_000, 91_000], "strand": ["+", "+", "+"]})
    gencode = pd.DataFrame({"chrom": "c1", "start0": [7_000, 51_000], "end1": [9_000, 51_500],
                            "strand": ["-", "-"], "gene_id": ["G1", "G2"], "gene_name": ["G1", "G2"],
                            "gene_type": "protein_coding", "transcript_id": ["t1", "t2"]})
    pos = mechanism.neighbour_position(loci, gencode)
    assert pos.to_dict() == {"anti": 1.0, "div": 1.0, "far": 0.0}


def test_specific_peptides():
    prot = pd.DataFrame({"peptide": ["A", "B", "C"], "lincORF_entries": ["e1", "e1,e2", "e2"],
                         "in_canonical_or_contaminant": ["False", "False", "True"]})
    hla = pd.DataFrame({"Peptide": ["H"], "lincORF_entries": ["e2"], "lincORF_specific": [True]})
    entries = pd.DataFrame({"entry": ["e1", "e2", "e2"], "orf_id": ["o1", "o2", "o3"]})   # o2, o3 identical
    hla_entries = pd.DataFrame({"entry": ["e2"], "orf_id": ["o2"]})
    p = mechanism.specific_peptides(prot, entries, hla, hla_entries)
    assert p.loc["o1", "n_specific_peptides"] == 2
    assert p.loc["o3", "n_specific_peptides"] == 1
    assert p.loc["o2", "n_hla_specific"] == 1
    assert p["copy_specific"].to_dict() == {"o1": True, "o2": True, "o3": False}


def test_class_table_phylocsf_weighting():
    feats = pd.DataFrame({c: [0.0, 1.0] for cols in mechanism.CLASSES.values() for c in cols})
    feats["m1_phylocsf"] = [0.0, 10.0]
    feats["copy_attribution"] = [1.0, 0.0]
    res = mechanism.class_table(feats.drop(columns=["m1_phylocsf_weighted"]), pd.Series([True, True]))
    assert res.loc[1, "pct_m1_phylocsf_weighted"] == pytest.approx(0.5)   # no attribution -> neutral
    assert res.loc[0, "pct_m1_phylocsf_weighted"] == pytest.approx(0.25)


def test_small_rna_hosts_dedupes_sources():
    loci = pd.DataFrame({"gene_id": ["host"], "chrom": "c1", "start0": [0], "end1": [1000], "strand": ["-"]})
    gencode = pd.DataFrame({"chrom": "c1", "start0": [100], "end1": [170], "strand": ["-"], "gene_id": ["S1"],
                            "gene_name": ["SNORD1"], "gene_type": ["snoRNA"], "transcript_id": ["t"]})
    wg = pd.DataFrame({"chrom": ["c1", "c1", "c1"], "start0": [100, 400, 500], "end1": [170, 470, 570],
                       "strand": ["-", "-", "+"]})
    sr = mechanism.small_rna_genes(gencode, wg)
    assert mechanism.small_rna_hosts(loci, sr).loc["host"] == 2.0
