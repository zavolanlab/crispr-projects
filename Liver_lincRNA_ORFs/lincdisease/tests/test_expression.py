import numpy as np
import pandas as pd
import pytest

from lincdisease import expression


def _regions():
    return pd.DataFrame({"feature_id": ["a", "a", "b"], "chrom": "chr1", "start0": [0, 100, 500],
                         "end1": [50, 150, 600]})


def _write_sums(d, regions, features):
    d.mkdir()
    regions.to_csv(d / "regions.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": ["s1", "s2"], **{f: [1.0, 2.0] for f in features}}).to_csv(
        d / "sums.0.tsv", sep="\t", index=False)


def test_read_sums_joins_added_features(tmp_path):
    r = _regions()
    _write_sums(tmp_path / "sums", r[r["feature_id"] == "a"], ["a"])
    _write_sums(tmp_path / "sums_b", r[r["feature_id"] == "b"], ["b"])
    sums = expression.read_sums(tmp_path, r)
    assert list(sums.columns) == ["sample_id", "a", "b"] and len(sums) == 2


def test_read_sums_raises_on_changed_regions(tmp_path):
    r = _regions()
    _write_sums(tmp_path / "sums", r, ["a", "b"])
    changed = r.copy()
    changed.loc[2, "end1"] = 700
    with pytest.raises(ValueError, match="changed"):
        expression.read_sums(tmp_path, changed)


def test_read_sums_raises_on_missing_feature(tmp_path):
    r = _regions()
    _write_sums(tmp_path / "sums", r[r["feature_id"] == "a"], ["a"])
    with pytest.raises(ValueError, match="no sums for 1 features"):
        expression.read_sums(tmp_path, r)


def test_read_sums_ignores_features_no_longer_quantified(tmp_path):
    r = _regions()
    _write_sums(tmp_path / "sums", r, ["a", "b"])
    sums = expression.read_sums(tmp_path, r[r["feature_id"] == "a"])
    assert list(sums.columns) == ["sample_id", "a"]


def test_size_factors_undo_composition():
    rng = np.random.default_rng(0)
    base = rng.lognormal(3, 1, size=(500, 1)) * np.ones((1, 20))
    depth = np.r_[np.ones(10), np.full(10, 2.0)]          # second half: every gene twice as high
    genes = pd.DataFrame(base * depth, columns=[f"s{i}" for i in range(20)])
    sf = expression.size_factors(genes)
    assert sf.iloc[10:].mean() / sf.iloc[:10].mean() == pytest.approx(2.0, rel=1e-6)
    assert np.exp(np.log(sf).mean()) == pytest.approx(1.0)


def test_read_sums_raises_on_samples_missing_from_added_features(tmp_path):
    r = _regions()
    _write_sums(tmp_path / "sums", r[r["feature_id"] == "a"], ["a"])
    d = tmp_path / "sums_b"
    d.mkdir()
    r[r["feature_id"] == "b"].to_csv(d / "regions.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": ["s1"], "b": [1.0]}).to_csv(d / "sums.0.tsv", sep="\t", index=False)   # no s2
    with pytest.raises(ValueError, match="1 samples x 1 features"):
        expression.read_sums(tmp_path, r)
