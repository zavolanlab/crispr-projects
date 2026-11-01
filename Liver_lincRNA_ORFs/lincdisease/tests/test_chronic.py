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


def _genes_qc(samples, seed=1):
    """Protein-coding gene sums with sample-specific depth and random technical QC."""
    rng = np.random.default_rng(seed)
    base = 2 ** rng.normal(5, 2, 300)
    depth = 2 ** rng.normal(0, 0.2, len(samples))
    genes = pd.DataFrame(np.outer(base, depth) * 2 ** rng.normal(0, 0.1, (300, len(samples))),
                         index=[f"G{i}" for i in range(300)], columns=samples)
    qc = pd.DataFrame({"star.number_of_input_reads_both": rng.integers(2e7, 4e7, len(samples)),
                       "star.uniquely_mapped_reads_%_both": rng.uniform(85, 95, len(samples))}, index=samples)
    return genes, qc


def test_d1_chronic_direction_and_combination():
    ncov, ph = _cohorts()
    genes, qc = _genes_qc(ncov.columns)
    res, sf = chronic.d1_chronic(ncov, ph, genes, qc, {}, genes.index)
    assert np.isclose(np.exp(np.log(sf).mean()), 1)
    assert res.loc["up", "D1_chronic"] > 0.8 and res.loc["down", "D1_chronic"] < -0.8
    assert abs(res.loc["flat", "D1_chronic"]) < 0.5
    # suppli has weight 0: reported, not combined
    assert res.loc["up", "chronic_fdr"] < 1e-3 and res.loc["up", "chronic_n_cohorts"] == 2
    assert res.loc["up", "hoang_r"] > 0.8 and res.loc["up", "suppli_r"] > 0.7
    assert not res["parent_driven"].any()


def test_combine_missing_cohort_renormalizes():
    cs = pd.DataFrame({"hoang_r": [0.4, np.nan], "hoang_p": [0.01, np.nan], "hoang_n": [50, np.nan],
                       "suppli_r": [0.2, 0.6], "suppli_p": [0.2, 0.01], "suppli_n": [40, 40]}, index=["a", "b"])
    out = chronic.combine(cs, weights={"hoang": 0.5, "suppli": 0.25})
    assert np.isclose(out.loc["a", "D1_chronic"], (0.5 * 0.4 + 0.25 * 0.2) / 0.75)
    assert np.isclose(out.loc["b", "D1_chronic"], 0.6)
    assert out.loc["b", "chronic_n_cohorts"] == 1


def test_combine_skips_zero_weight_cohorts():
    cs = pd.DataFrame({"hoang_r": [0.4], "hoang_p": [0.01], "hoang_n": [50],
                       "suppli_r": [0.9], "suppli_p": [1e-6], "suppli_n": [40]}, index=["a"])
    out = chronic.combine(cs)          # default WEIGHTS: suppli 0
    assert np.isclose(out.loc["a", "D1_chronic"], 0.4) and out.loc["a", "chronic_n_cohorts"] == 1
