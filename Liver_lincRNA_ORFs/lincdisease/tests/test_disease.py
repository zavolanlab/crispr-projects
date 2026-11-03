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
    sf = pd.Series(np.exp(rng.normal(0, 0.2, len(tum))), index=tum["sample_id"])
    res = disease.d2_survival(expr, s, cdr, sf, ncov=ncov)
    assert res.loc["risk", "HR_OS"] > 2
    assert res.loc["risk", "D2"] > 5
    assert abs(res.loc["null", "D2"]) < 3


def test_disease_coverage_counts_missing_subcomponents():
    from lincdisease import integrate
    rng = np.random.default_rng(1)
    stats = pd.DataFrame(rng.normal(size=(50, 9)), columns=["D1_HCC", "D1_chronic", "D2", "D3", "D4", "D5",
                                                           "gtex_liver_median", "liver_specificity",
                                                           "cell_line_log2TPM"])
    stats.loc[0, "D1_chronic"] = np.nan                     # 0.3 of D1 (weight 0.25) missing
    res = integrate.disease_score(stats, pd.Series(True, index=stats.index))
    assert res.loc[0, "disease_coverage"] == pytest.approx(1 - 0.25 * 0.3)
    assert res.loc[1, "disease_coverage"] == pytest.approx(1.0)


def test_d1_hcc_is_the_paired_statistic():
    s = _samples()
    rng = np.random.default_rng(3)
    expr = pd.DataFrame(rng.normal(size=(1, len(s))), index=["x"], columns=s["sample_id"])
    expr.loc["x", s.loc[s["source"] == "gtex", "sample_id"]] += 5      # GTEx batch: no effect on D1 HCC
    res = disease.d1_hcc(expr, s)
    assert res.loc["x", "D1_HCC"] == res.loc["x", "r_paired"]
    assert res.loc["x", "r_vs_gtex"] < -0.9


def test_hcc_and_chronic_scores_and_direction():
    from lincdisease import integrate
    rng = np.random.default_rng(2)
    cols = ["D1_HCC", "D1_chronic", "D2", "D3", "D4", "D4_hcc", "D4_chronic", "D5", "D5_hcc", "D5_chronic",
            "gtex_liver_median", "liver_specificity", "cell_line_log2TPM"]
    stats = pd.DataFrame(rng.normal(size=(200, len(cols))), columns=cols)
    stats.loc[0, ["D1_HCC", "D2", "D3"]] = [10, 10, 10]          # strong and up in every HCC component
    stats.loc[1, ["D1_HCC", "D2", "D3"]] = [-10, -10, -10]
    stats.loc[2, "D1_chronic"] = -10                               # strong, down with disease
    res = integrate.disease_score(stats, pd.Series(True, index=stats.index))
    assert res.loc[0, "hcc_direction"] == "oncogene-like" and res.loc[1, "hcc_direction"] == "suppressor-like"
    # magnitude is ranked: down-regulated loci score as high as up-regulated ones
    assert res.loc[1, "comp_D1_HCC"] == pytest.approx(res.loc[0, "comp_D1_HCC"])
    assert res.loc[2, "chronic_direction"] == "down" and res.loc[2, "comp_D1_chronic"] > 0.99
    assert res["hcc_score"].between(0, 1).all() and res["chronic_score"].between(0, 1).all()
