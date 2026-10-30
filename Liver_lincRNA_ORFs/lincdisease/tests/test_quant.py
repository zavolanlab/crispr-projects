import numpy as np
import pandas as pd
import pytest

from lincdisease import quant


def test_merge_intervals():
    df = pd.DataFrame({"feature_id": ["a"] * 3 + ["b"], "chrom": ["c1"] * 4,
                       "start0": [10, 15, 40, 0], "end1": [20, 30, 50, 5]})
    m = quant.merge_intervals(df)
    assert m[m.feature_id == "a"][["start0", "end1"]].values.tolist() == [[10, 30], [40, 50]]


def test_subtract_intervals_splits():
    a = pd.DataFrame({"feature_id": ["a"], "chrom": ["c1"], "start0": [0], "end1": [100]})
    b = pd.DataFrame({"chrom": ["c1", "c1", "c2"], "start0": [10, 50, 0], "end1": [20, 200, 100]})
    r = quant.subtract_intervals(a, b)
    assert r[["start0", "end1"]].values.tolist() == [[0, 10], [20, 50]]


def _exons(rows):
    return pd.DataFrame(rows, columns=["chrom", "start0", "end1", "strand", "gene_id", "gene_name",
                                       "gene_type", "transcript_id"])


def test_build_regions_keeps_own_gene_and_removes_antisense():
    linc = pd.DataFrame({"gene_id": ["L1", "L1"], "chrom": ["c1", "c1"], "start0": [0, 300], "end1": [200, 400]})
    gencode = _exons([
        ("c1", 0, 200, "+", "ENSG1.1", "LINCX", "lncRNA", "T1"),      # the locus itself
        ("c1", 150, 250, "-", "ENSG2.3", "PCG", "protein_coding", "T2"),  # antisense coding exon
        ("c1", 1000, 1300, "+", "ENSG3.1", "CTRL", "lncRNA", "T3"),
    ])
    loci = pd.DataFrame({"name_source_gene_ids": ["ENSG1.1"], "lncRNA_gene_ids": ["ENSG1.1"]},
                        index=pd.Index(["L1"], name="locus_id"))
    regions, feats = quant.build_regions(linc, gencode, loci, ["CTRL"], min_quant_bases=100)
    r = regions[regions.feature_id == "L1"][["start0", "end1"]].values.tolist()
    assert r == [[0, 150], [300, 400]]
    assert feats.loc["L1", "exonic_bases"] == 300
    assert feats.loc["L1", "quant_bases"] == 250
    assert feats.loc["ctrl:CTRL", "quant_bases"] == 300
    assert feats["quantifiable"].all()


def test_build_regions_unknown_control():
    linc = pd.DataFrame({"gene_id": ["L1"], "chrom": ["c1"], "start0": [0], "end1": [200]})
    loci = pd.DataFrame({"name_source_gene_ids": [""], "lncRNA_gene_ids": [""]}, index=["L1"])
    with pytest.raises(ValueError):
        quant.build_regions(linc, _exons([]), loci, ["NOPE"])


class FakeBW:
    def chroms(self):
        return {"c1": 1000}

    def stats(self, c, s, e, type, exact):
        return [float(e - s)]          # coverage 1 everywhere


def test_sum_regions_and_normalize():
    regions = pd.DataFrame({"feature_id": ["a", "a", "b"], "chrom": ["c1", "c1", "c9"],
                            "start0": [0, 100, 0], "end1": [50, 150, 10]})
    s = quant.sum_regions(FakeBW(), regions)
    assert s.to_dict() == {"a": 100.0, "b": 0.0}
    sums = pd.DataFrame({"sample_id": ["s1"], "a": [100.0], "b": [0.0]})
    feats = pd.DataFrame({"quant_bases": [100, 10]}, index=["a", "b"])
    samples = pd.DataFrame({"sample_id": ["s1"], "auc": [1e9]})
    n = quant.normalize(sums, feats, samples)
    assert n.loc["a", "s1"] == pytest.approx(1.0)
    assert n.loc["b", "s1"] == 0.0
