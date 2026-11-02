import numpy as np
import pandas as pd
import pytest

from lincdisease import scoring


def test_percentile_midrank_and_missing():
    bg = [1, 2, 3, 4]
    p = scoring.percentile_vs_background([0, 2, 2.5, 5, np.nan], bg)
    assert p[0] == 0.0
    assert p[1] == pytest.approx(1.5 / 4)      # one below, one tie
    assert p[2] == pytest.approx(2 / 4)
    assert p[3] == 1.0
    assert np.isnan(p[4])


def test_typical_value_is_half():
    bg = np.arange(1001)
    assert scoring.percentile_vs_background([500], bg)[0] == pytest.approx(0.5)


def test_weighted_percentile_counts_each_locus_once():
    # locus A has three ORFs with value 1, locus B one ORF with value 3: weighted, 1 is half the background
    bg = [1, 1, 1, 3]
    w = scoring.locus_weights(["A", "A", "A", "B"])
    assert w.tolist() == pytest.approx([1 / 3, 1 / 3, 1 / 3, 1])
    assert scoring.percentile_vs_background([1, 2], bg, w).tolist() == pytest.approx([0.25, 0.5])
    assert scoring.percentile_vs_background([1, 2], bg).tolist() == pytest.approx([0.375, 0.75])


def test_signed_percentile_ranks_magnitude():
    pct, sign = scoring.signed_percentile_vs_background([-3, 3, 0.1], [-1, 1, 2, -2])
    assert pct[0] == pct[1] == 1.0
    assert list(sign) == [-1, 1, 1]


def test_composite_shrinks_missing_to_half():
    p = pd.DataFrame({"a": [1.0, 1.0, np.nan], "b": [1.0, np.nan, np.nan]})
    res = scoring.composite(p, {"a": 0.75, "b": 0.25})
    assert res["score"].tolist() == pytest.approx([1.0, 0.75 + 0.25 * 0.5, 0.5])
    assert res["coverage"].tolist() == pytest.approx([1.0, 0.75, 0.0])


def test_composite_rejects_bad_weights():
    with pytest.raises(ValueError):
        scoring.composite(pd.DataFrame({"a": [0.1]}), {"a": -1})


def test_class_scores_standardized_and_dominance():
    feats = pd.DataFrame({"f1": [0.9, 0.5], "f2": [0.9, 0.5], "g1": [0.3, 0.55]})
    scores, cov = scoring.class_scores(feats, {"A": ["f1", "f2"], "B": ["g1", "missing"]})
    assert scores.loc[0, "A"] == pytest.approx(0.9)
    assert scores.loc[0, "B"] == pytest.approx(0.4)      # missing feature -> 0.5
    assert cov.loc[0, "B"] == pytest.approx(0.5)
    z = scoring.standardized(scores, np.array([True, True]))
    assert z.loc[0, "A"] > 0 > z.loc[1, "A"]
    assert z.loc[0, "A"] == pytest.approx(-z.loc[1, "A"])


def test_standardized_is_rank_based_and_skips_classes_without_data():
    # class R: one rare binary feature (1 of 20 rows); class C: continuous. Rank-based standardization keeps the rare
    # class within the range of the continuous one
    scores = pd.DataFrame({"R": [1.0] + [0.47] * 19, "C": np.linspace(0, 1, 20)})
    cov = pd.DataFrame({"R": [1.0] * 20, "C": [1.0] * 19 + [0.0]})
    z = scoring.standardized(scores, np.ones(20, bool), cov)
    assert z.loc[0, "R"] < 2.0                             # (score - mean) / SD would give 4.2
    assert np.isnan(z.loc[19, "C"])                       # no data for C: no standardized score
    assert z["C"].iloc[:19].median() == pytest.approx(0.0, abs=1e-9)
    z = pd.DataFrame({"A": [2.0, 0.6, -0.5], "B": [0.5, 0.3, -1.0]})
    # ahead by 1.5 SD; ahead by only 0.3; ahead by 0.5 but below the mean
    assert scoring.dominant_class(z, margin=0.5).tolist() == ["A", "none", "none"]


def test_consensus_direction():
    pct = pd.DataFrame({"D1": [0.9, 0.9, 0.9, 0.2], "D2": [0.95, 0.95, 0.1, 0.9]})
    sign = pd.DataFrame({"D1": [1, 1, 1, 1], "D2": [1, -1, -1, np.nan]})
    d = scoring.consensus_direction(pct, sign, 0.75, labels=("oncogene-like", "suppressor-like"))
    assert d.tolist() == ["oncogene-like", "inconsistent", "oncogene-like", "none"]
