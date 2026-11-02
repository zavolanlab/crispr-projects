import numpy as np
import pandas as pd
import pytest

from lincdisease import hcc


def test_sample_labels():
    assert hcc.sample_labels("GSE77314", "GSM1: Patient 12 Tumor; Homo sapiens", {}) == ("12", "tumour")
    assert hcc.sample_labels("GSE94660", "GSM2: Non-neoplastic liver tissue of P105; Homo", {}) == ("P105", "normal")
    assert hcc.sample_labels("GSE56545", "", {"source_name": "HCC-C11 Non-recurrent Tumor's Adjacent Normal"}) == \
        ("HCC-C11", "normal")
    assert hcc.sample_labels("GSE77509", "", {"patient id": "#3", "tissue": "PVTT"}) == ("#3", None)


def _cohort(cohort, n, up):
    """n patients; feature 'x' up 4-fold in tumours if `up`, 'flat' unchanged; GSE77509-style two runs per sample."""
    rows, cols = [], {}
    rng = np.random.default_rng(len(cohort))
    for i in range(n):
        for tissue in ("normal", "tumour"):
            exp = f"{cohort}_{i}_{tissue}"
            for run in range(2):
                sid = f"{exp}_r{run}"
                rows.append(dict(sample_id=sid, cohort=cohort, experiment_acc=exp, patient=f"{cohort}:{i}",
                                 tissue=tissue))
                base = rng.uniform(1, 2)
                cols[sid] = [base * (4 if up and tissue == "tumour" else 1), base]
    return pd.DataFrame(rows), pd.DataFrame(cols, index=["x", "flat"])


def test_replication_combines_cohorts():
    s1, n1 = _cohort("GSE77314", 12, True)
    s2, n2 = _cohort("GSE94660", 8, True)
    res = hcc.replication(pd.concat([n1, n2], axis=1), pd.concat([s1, s2], ignore_index=True))
    assert res.loc["x", "GSE77314_r"] == pytest.approx(1.0) and res.loc["x", "GSE77314_pairs"] == 12
    assert res.loc["x", "hcc_rep_r"] == pytest.approx(1.0)
    assert res.loc["x", "hcc_rep_pairs"] == 20 and res.loc["x", "hcc_rep_n_cohorts"] == 2
    assert res.loc["x", "hcc_rep_p"] < 1e-3
    assert abs(res.loc["flat", "hcc_rep_r"]) < 0.6
