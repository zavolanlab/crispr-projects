import numpy as np
import pandas as pd
import pytest

from lincdisease import nmd


def _data():
    cond = ["control"] * 4 + ["UPF1_kd"] * 2 + ["SMG6_SMG7_kd"] * 2
    sid = [f"s{i}" for i in range(len(cond))]
    samples = pd.DataFrame({"sample_id": sid, "condition": cond})
    ncov = pd.DataFrame([[1.0] * 4 + [4.0] * 2 + [4.0] * 2,          # NMD target: up 4-fold in both
                         [2.0] * 8,                                   # unchanged
                         [0.0] * 7 + [1.0]],                          # barely detected
                        index=["target", "flat", "rare"], columns=sid)
    rng = np.random.default_rng(0)
    genes = pd.DataFrame(rng.uniform(50, 100, size=(200, 1)) * np.ones((1, len(sid))), columns=sid)
    return ncov, samples, genes


def test_sensitivity_direction_and_detection():
    ncov, samples, genes = _data()
    res = nmd.sensitivity(ncov, samples, genes, pseudocount=0.0001)
    assert res.loc["target", "nmd_log2fc"] == pytest.approx(2.0, abs=1e-3)
    assert res.loc["target", "log2fc_UPF1_kd"] == pytest.approx(res.loc["target", "log2fc_SMG6_SMG7_kd"])
    assert res.loc["flat", "nmd_log2fc"] == pytest.approx(0.0, abs=1e-6)
    assert np.isnan(res.loc["rare", "nmd_log2fc"]) and res.loc["rare", "detected_fraction"] == pytest.approx(1 / 8)


def test_size_factors_remove_sample_depth():
    ncov, samples, genes = _data()
    # one control sample sequenced twice as deep: features and genes double, the result must not change
    ncov2, genes2 = ncov.copy(), genes.copy()
    ncov2["s0"] *= 2
    genes2["s0"] *= 2
    a = nmd.sensitivity(ncov, samples, genes, pseudocount=0.0001)
    b = nmd.sensitivity(ncov2, samples, genes2, pseudocount=0.0001)
    assert b.loc["target", "nmd_log2fc"] == pytest.approx(a.loc["target", "nmd_log2fc"], abs=1e-3)
