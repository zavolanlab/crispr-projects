import numpy as np
import pandas as pd
import pytest

from lincdisease import disease


def test_rank_biserial_extremes():
    r, p = disease.rank_biserial([5, 6, 7], [1, 2, 3])
    assert r == pytest.approx(1.0)
    r, _ = disease.rank_biserial([1, 2, 3], [5, 6, 7])
    assert r == pytest.approx(-1.0)


def test_paired_rank_biserial():
    r, _ = disease.paired_rank_biserial([1, 2, 3, 0])
    assert r == pytest.approx(1.0)
    r, _ = disease.paired_rank_biserial([1, -1])
    assert r == pytest.approx(0.0)


def _samples(n_pat=30, n_gtex=30):
    rows = []
    for i in range(n_pat):
        bc = f"TCGA-AA-{i:04d}"
        rows.append(dict(source="tcga", project="LIHC", sample_id=f"t{i}", tcga_barcode=bc + "-01A",
                         cgc_sample_sample_type="Primary Tumor"))
        if i < 10:
            rows.append(dict(source="tcga", project="LIHC", sample_id=f"n{i}", tcga_barcode=bc + "-11A",
                             cgc_sample_sample_type="Solid Tissue Normal"))
    for j in range(n_gtex):
        rows.append(dict(source="gtex", project="LIVER", sample_id=f"g{j}"))
    return pd.DataFrame(rows)


def test_d1_detects_up_and_down():
    s = _samples()
    rng = np.random.default_rng(0)
    cols = s["sample_id"]
    tumour = cols.str.startswith("t").to_numpy()
    up = np.where(tumour, 5.0, 1.0) + rng.normal(0, 0.1, len(cols))
    flat = rng.normal(3, 1, len(cols))
    expr = pd.DataFrame([up, -up, flat], index=["up", "down", "flat"], columns=cols)
    res = disease.d1_hcc(expr, s)
    assert res.loc["up", "D1_HCC"] == pytest.approx(1.0)
    assert res.loc["down", "D1_HCC"] == pytest.approx(-1.0)
    assert abs(res.loc["flat", "D1_HCC"]) < 0.5
    assert res.loc["up", "fdr_paired"] < 0.05
    assert res.loc["up", "n_pairs"] == 10


def test_d2_hazard_direction():
    s = _samples(n_pat=200, n_gtex=0)
    rng = np.random.default_rng(1)
    tum = s[s["cgc_sample_sample_type"] == "Primary Tumor"].reset_index(drop=True)
    x = rng.normal(0, 1, len(tum))
    t_event = rng.exponential(1 / np.exp(1.0 * x)) * 1000
    cdr = pd.DataFrame({"bcr_patient_barcode": tum["tcga_barcode"].str[:12], "type": "LIHC",
                        "ajcc_pathologic_tumor_stage": "Stage I", "age_at_initial_pathologic_diagnosis": 60,
                        "OS": 1, "OS.time": t_event, "PFI": 1, "PFI.time": t_event})
    cdr.loc[::2, "ajcc_pathologic_tumor_stage"] = "Stage IIIA"
    cdr["age_at_initial_pathologic_diagnosis"] = rng.integers(40, 80, len(cdr))
    expr = pd.DataFrame([x, rng.normal(0, 1, len(tum))], index=["risk", "null"], columns=tum["sample_id"])
    ncov = expr * 0 + 1
    res = disease.d2_survival(expr, s, cdr, ncov=ncov)
    assert res.loc["risk", "HR_OS"] > 2
    assert res.loc["risk", "D2"] > 5
    assert abs(res.loc["null", "D2"]) < 3
