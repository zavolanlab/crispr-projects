import numpy as np
import pandas as pd

from lincdisease import chronic


def _cohorts(seed=0):
    rng = np.random.default_rng(seed)
    man = pd.DataFrame({
        "project": ["SRP197353"] * 40 + ["SRP186450"] * 20 + ["SRP078450"] * 12,
        "sample_attributes": [f"fibrosis stage;;{i % 5}|nafld activity score;;{i % 7}" for i in range(40)]
        + [f"disease;;{d}" for d in ["healthy", "obese", "NAFLD", "NASH"] * 5]
        + [f"infection status;;{s}|treatment;;None" for s in ["HCV", "non-HCV"] * 5]
        + ["infection status;;HCV|treatment;;pegIFNa/ribavirin 4h"] * 2,
    }, index=[f"s{i}" for i in range(72)])
    ph = chronic.phenotypes(man)
    sev = ph["fibrosis_stage"].fillna(ph["group"].map({"case": 3, "control": 0}))
    up = 2 ** (1 + sev.fillna(0) + rng.normal(0, 0.3, len(sev)))
    flat = pd.Series(2 ** (1 + rng.normal(0, 0.3, len(sev))), index=sev.index)
    ncov = pd.DataFrame([up, 1 / up, flat], index=["up", "down", "flat"])
    return ncov, ph


def test_phenotypes():
    _, ph = _cohorts()
    assert ph["cohort"].value_counts().to_dict() == {"hoang": 40, "suppli": 20, "hcv": 12}
    assert ph.loc["s40":"s59", "group"].value_counts().to_dict() == {"control": 10, "case": 10}
    assert ph.loc[["s70", "s71"], "group"].isna().all()       # IFN-treated: not used
    assert ph.loc["s3", "fibrosis_stage"] == 3 and ph.loc["s3", "nas"] == 3


def test_d1_chronic_direction_and_combination():
    ncov, ph = _cohorts()
    res = chronic.d1_chronic(ncov, ph)
    assert res.loc["up", "D1_chronic"] > 0.8 and res.loc["down", "D1_chronic"] < -0.8
    assert abs(res.loc["flat", "D1_chronic"]) < 0.5
    assert res.loc["up", "chronic_fdr"] < 1e-3 and res.loc["up", "chronic_n_cohorts"] == 3
    assert res.loc["up", "hoang_r"] > 0.9 and res.loc["up", "suppli_r"] == 1.0


def test_combine_missing_cohort_renormalizes():
    cs = pd.DataFrame({"hoang_r": [0.4, np.nan], "hoang_p": [0.01, np.nan], "hoang_n": [50, np.nan],
                       "suppli_r": [0.2, 0.6], "suppli_p": [0.2, 0.01], "suppli_n": [40, 40]}, index=["a", "b"])
    out = chronic.combine(cs)
    assert np.isclose(out.loc["a", "D1_chronic"], (0.5 * 0.4 + 0.25 * 0.2) / 0.75)
    assert np.isclose(out.loc["b", "D1_chronic"], 0.6)
    assert out.loc["b", "chronic_n_cohorts"] == 1
