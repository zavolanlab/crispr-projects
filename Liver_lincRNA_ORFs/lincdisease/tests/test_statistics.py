"""Shared statistics: cohort combination (disease.stouffer / weighted_mean), partial rank correlation
(chronic.partial_rank), D6 / tumour specificity (integrate), the D4 window (genetics.d4)."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from lincdisease import chronic, disease, genetics, integrate


def test_stouffer_single_cohort_returns_its_own_z():
    r, p, n = np.array([[0.4, 0.0]]), np.array([[0.01, 0.5]]), np.array([[30, 30]])
    z, pc = disease.stouffer(r, p, n, np.array([[True, False]]))
    assert z[0] == pytest.approx(stats.norm.isf(0.005)) and pc[0] == pytest.approx(0.01)


def test_stouffer_weights_sqrt_n_and_skips_cells():
    r = np.array([[1.0, -1.0], [1.0, np.nan], [np.nan, np.nan]])
    p = np.array([[0.05, 0.05], [0.05, np.nan], [np.nan, np.nan]])
    n = np.array([[100, 25], [100, 25], [100, 25]])
    ok = np.isfinite(r)
    z, pc = disease.stouffer(r, p, n, ok)
    z1 = stats.norm.isf(0.025)
    assert z[0] == pytest.approx((10 * z1 - 5 * z1) / np.sqrt(125))
    assert z[1] == pytest.approx(z1)            # the NaN cohort does not count
    assert np.isnan(z[2]) and np.isnan(pc[2])


def test_weighted_mean():
    r = np.array([[0.2, 0.6], [0.2, np.nan], [np.nan, np.nan]])
    m = disease.weighted_mean(r, np.array([1.0, 3.0]), np.isfinite(r))
    assert m[0] == pytest.approx(0.5) and m[1] == pytest.approx(0.2) and np.isnan(m[2])


def test_partial_rank_equals_regression_partial_correlation():
    rng = np.random.default_rng(3)
    n = 60
    c = pd.DataFrame({"c": rng.normal(size=n)}, index=[f"s{i}" for i in range(n)])
    y = pd.Series(c["c"].to_numpy() + rng.normal(size=n), index=c.index)
    expr = pd.DataFrame([c["c"].to_numpy() + rng.normal(size=n), rng.normal(size=n)], index=["f1", "f2"],
                        columns=c.index)
    out = chronic.partial_rank(expr, y, c)
    rc = c["c"].rank().to_numpy()
    X = np.column_stack([np.ones(n), rc])

    def resid(v):
        return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
    for f in expr.index:
        r = np.corrcoef(resid(expr.loc[f].rank().to_numpy()), resid(y.rank().to_numpy()))[0, 1]
        t = r * np.sqrt((n - 3) / (1 - r ** 2))
        assert out.loc[f, "stat"] == pytest.approx(r)
        assert out.loc[f, "p"] == pytest.approx(2 * stats.t.sf(abs(t), n - 3))


def _samples():
    rows = [("tcga", "LIHC", "Primary Tumor", "t1"), ("tcga", "LIHC", "Primary Tumor", "t2"),
            ("tcga", "LIHC", "Solid Tissue Normal", "n1"), ("gtex", "LIVER", None, "g1"),
            ("gtex", "LIVER", None, "g2"), ("gtex", "BRAIN", None, "g3"), ("gtex", "TESTIS", None, "g4")]
    return pd.DataFrame(rows, columns=["source", "project", "cgc_sample_sample_type", "sample_id"])


def test_liver_context_and_tumour_specificity():
    expr = pd.DataFrame({"t1": [5.0], "t2": [7.0], "n1": [0.0], "g1": [2.0], "g2": [4.0], "g3": [1.0],
                         "g4": [9.0]}, index=["f"])
    d6 = integrate.liver_context(expr, _samples(), pd.Series({"f": 1.5}))
    assert d6.loc["f", "gtex_liver_median"] == 3.0
    assert d6.loc["f", "liver_specificity"] == 3.0 - 9.0         # testis is a GTEx tissue for D6
    assert d6.loc["f", "cell_line_log2TPM"] == 1.5
    # tumour median 6 (adjacent normals left out) minus the largest GTEx median without testis (liver 3)
    assert integrate.tumour_specificity(expr, _samples())["f"] == 3.0


def test_d4_window_edges():
    spans = pd.DataFrame({"chrom": ["chr1"], "start0": [100_000], "end1": [100_500]}, index=["L"])

    def d4_at(pos):
        assoc = pd.DataFrame({"chrom": ["chr1"], "pos": [pos], "mlog10p": [10.0], "trait": ["cirrhosis of liver"],
                              "trait_class": ["disease"], "snp": ["rs1"]})
        return genetics.d4(spans, assoc, flank=1000).loc["L", "D4"]
    # the span is the 1-based positions 100001..100500; +- 1 kb: 99001..101500
    assert d4_at(99_001) == 10.0 and d4_at(101_500) == 10.0
    assert d4_at(99_000) == 0.0 and d4_at(101_501) == 0.0
