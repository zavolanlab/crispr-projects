"""
Disease components (LOGIC.md, Step 3) from the recount3 quantification (quant.normalize: features x samples).

Expression is log2(ncov + pseudocount). Effect sizes are rank-based, so they are comparable between loci with
very different expression levels and robust to the many zeros of lowly expressed lincRNAs:

D1 HCC   rank-biserial correlation r (-1..1) of
           paired   TCGA-LIHC primary tumour vs adjacent normal of the same patient (matched-pairs, from the
                    Wilcoxon signed-rank test)
           GTEx     TCGA-LIHC primary tumours vs GTEx liver (Mann-Whitney; AUC = (r + 1) / 2)
         statistic D1_HCC = mean of the two r (signed: > 0 up in tumour); p-values BH-adjusted over features.
D2       Cox proportional hazards of overall survival (OS) and progression-free interval (PFI) (TCGA-CDR
         endpoints) on standardized expression, adjusted for age and stage (I-II vs III-IV), one primary tumour
         per patient; features expressed (ncov > 0) in < `min_expressed` of tumours are not tested.
         statistic D2 = mean of the Wald z of OS and PFI (signed: > 0 = high expression, worse outcome).
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def log_expr(ncov, pseudocount):
    return np.log2(ncov + pseudocount)


def bh(p):
    p = pd.Series(p, dtype=float)
    out = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def rank_biserial(x, y):
    """Mann-Whitney rank-biserial correlation of x vs y (> 0: x larger) and two-sided p."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u / (len(x) * len(y)) - 1, p


def paired_rank_biserial(diff):
    """Matched-pairs rank-biserial correlation of paired differences (zeros dropped) and Wilcoxon p."""
    d = np.asarray(diff, float)
    d = d[d != 0]
    if len(d) == 0:
        return 0.0, 1.0
    ranks = stats.rankdata(np.abs(d))
    r = (ranks[d > 0].sum() - ranks[d < 0].sum()) / ranks.sum()
    return r, stats.wilcoxon(d).pvalue


def tcga_groups(samples):
    """Primary tumours (one per patient), adjacent normals, and GTEx liver sample ids, plus patient barcodes."""
    t = samples[samples["source"] == "tcga"].copy()
    t["patient"] = t["tcga_barcode"].str[:12]
    tum = t[t["cgc_sample_sample_type"] == "Primary Tumor"].sort_values("tcga_barcode").drop_duplicates("patient")
    nor = t[t["cgc_sample_sample_type"] == "Solid Tissue Normal"].sort_values("tcga_barcode").drop_duplicates("patient")
    gtex_liver = samples[(samples["source"] == "gtex") & (samples["project"] == "LIVER")]
    return tum, nor, gtex_liver


def d1_hcc(expr, samples):
    tum, nor, liver = tcga_groups(samples)
    pairs = tum.merge(nor, on="patient", suffixes=("_t", "_n"))
    et = expr[tum["sample_id"]]
    eg = expr[liver["sample_id"]]
    diff = expr[pairs["sample_id_t"]].to_numpy() - expr[pairs["sample_id_n"]].to_numpy()
    rows = []
    for i, fid in enumerate(expr.index):
        r_p, p_p = paired_rank_biserial(diff[i])
        r_g, p_g = rank_biserial(et.iloc[i], eg.iloc[i])
        rows.append((fid, diff[i].mean(), r_p, p_p, et.iloc[i].mean() - eg.iloc[i].mean(), (r_g + 1) / 2, r_g, p_g))
    out = pd.DataFrame(rows, columns=["feature_id", "log2FC_paired", "r_paired", "p_paired",
                                      "log2FC_vs_gtex", "auc_vs_gtex", "r_vs_gtex", "p_vs_gtex"]).set_index("feature_id")
    out["fdr_paired"] = bh(out["p_paired"])
    out["fdr_vs_gtex"] = bh(out["p_vs_gtex"])
    out["n_pairs"] = len(pairs)
    out["n_tumour"] = len(tum)
    out["n_gtex_liver"] = len(liver)
    out["D1_HCC"] = out[["r_paired", "r_vs_gtex"]].mean(axis=1)
    return out


STAGE_LATE = ("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV", "Stage IVA", "Stage IVB")
STAGE_EARLY = ("Stage I", "Stage II")


def survival_table(samples, cdr):
    """Primary tumours with OS/PFI endpoints, age and late-stage indicator (patients without stage dropped)."""
    tum, _, _ = tcga_groups(samples)
    c = cdr[cdr["type"] == "LIHC"].copy()
    c = c[c["ajcc_pathologic_tumor_stage"].isin(STAGE_LATE + STAGE_EARLY)]
    c["late_stage"] = c["ajcc_pathologic_tumor_stage"].isin(STAGE_LATE).astype(float)
    c["age"] = pd.to_numeric(c["age_at_initial_pathologic_diagnosis"], errors="coerce")
    for col in ("OS", "OS.time", "PFI", "PFI.time"):
        c[col] = pd.to_numeric(c[col], errors="coerce")
    t = tum.merge(c, left_on="patient", right_on="bcr_patient_barcode")
    return t.dropna(subset=["age"])


def _cox(df, endpoint):
    from lifelines import CoxPHFitter
    d = df[["x", "age", "late_stage", endpoint, f"{endpoint}.time"]].dropna()
    d = d[d[f"{endpoint}.time"] > 0]
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(d, duration_col=f"{endpoint}.time", event_col=endpoint)
    s = cph.summary.loc["x"]
    return s["exp(coef)"], s["z"], s["p"], len(d)


def d2_survival(expr, samples, cdr, min_expressed=0.2, ncov=None):
    """Cox models per feature; `ncov` (linear scale) decides which features are expressed enough to test."""
    tab = survival_table(samples, cdr)
    e = expr[tab["sample_id"]]
    expressed = (ncov[tab["sample_id"]] > 0).mean(axis=1) if ncov is not None else pd.Series(1.0, index=e.index)
    rows = []
    for fid in e.index:
        rec = {"feature_id": fid, "frac_tumours_expressed": expressed[fid]}
        x = e.loc[fid].to_numpy()
        if expressed[fid] >= min_expressed and np.std(x) > 0:
            df = tab[["age", "late_stage", "OS", "OS.time", "PFI", "PFI.time"]].assign(x=(x - x.mean()) / x.std())
            for ep in ("OS", "PFI"):
                try:
                    rec[f"HR_{ep}"], rec[f"z_{ep}"], rec[f"p_{ep}"], rec[f"n_{ep}"] = _cox(df, ep)
                except Exception:          # non-convergence: leave missing
                    pass
        rows.append(rec)
    out = pd.DataFrame(rows).set_index("feature_id")
    for ep in ("OS", "PFI"):
        for col in (f"HR_{ep}", f"z_{ep}", f"p_{ep}", f"n_{ep}"):
            if col not in out:
                out[col] = np.nan
        out[f"fdr_{ep}"] = bh(out[f"p_{ep}"])
    out["D2"] = out[["z_OS", "z_PFI"]].mean(axis=1)
    return out
