"""
Disease components (LOGIC.md, Step 3) from the recount3 quantification of TCGA-LIHC and GTEx (features x samples,
ncov scaled by the size factors of expression.py).

Expression is log2(ncov + pseudocount). Effect sizes are rank-based, so they are comparable between loci with
very different expression levels; zeros (pairs without a difference) are kept in the ranking, so a locus that
is undetected in most samples cannot reach an extreme r from its few detected ones:

D1 HCC   matched-pairs rank-biserial correlation r (-1..1) of TCGA-LIHC primary tumour vs adjacent normal of the
         same patient, from the Wilcoxon signed-rank test with zero differences kept (Pratt; signed: > 0 up in
         tumour); features detected (ncov > 0) in < MIN_DETECTED of the paired samples are not tested; p-values
         BH-adjusted over features. Reported, not scored: tumours vs GTEx liver (Mann-Whitney r; AUC = (r + 1) / 2), as TCGA adjacent
         normals and GTEx liver, both non-tumour liver, already differ for most loci (study batch).
D2       Cox proportional hazards of overall survival (OS) and progression-free interval (PFI) (TCGA-CDR
         endpoints) on standardized expression, adjusted for age, stage (I-II vs III-IV) and the log size factor
         of the tumour (the loss of hepatocyte transcripts in less differentiated tumours, which itself predicts
         survival), one primary tumour per patient; features expressed (ncov > 0) in < `min_expressed` of tumours
         are not tested. statistic D2 = mean of the Wald z of OS and PFI (signed: > 0 = high expression, worse
         outcome).
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


MIN_DETECTED = 0.2      # features detected (ncov > 0) in fewer of the samples are not tested (D1 HCC, D2)


def log_expr(ncov, pseudocount):
    """log2(ncov + pseudocount), the expression scale of every statistic."""
    return np.log2(ncov + pseudocount)


def bh(p):
    """Benjamini-Hochberg FDR of the non-missing p-values (missing stay missing)."""
    p = pd.Series(p, dtype=float)
    out = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def weighted_mean(r, w, ok):
    """Row-wise mean of r (features x cohorts) weighted by w over the cells where `ok`; NaN where no weight."""
    w = np.where(ok, w, 0.0)
    wsum = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(wsum > 0, (np.where(ok, r, 0.0) * w).sum(axis=1) / wsum, np.nan)


def stouffer(r, p, n, ok):
    """
    Row-wise Stouffer combination over the cells where `ok` (features x cohorts): z of each two-sided p with the
    sign of its statistic r, weights sqrt(n). Returns the combined z and its two-sided p (NaN where no cell is ok).
    """
    z = np.sign(r) * stats.norm.isf(np.clip(p, 1e-300, 1) / 2)
    sw = np.where(ok, np.sqrt(n), 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        zc = np.nansum(np.where(ok, z * sw, 0), axis=1) / np.sqrt((sw ** 2).sum(axis=1))
    zc = np.where(ok.any(axis=1), zc, np.nan)
    return zc, np.where(np.isfinite(zc), 2 * stats.norm.sf(np.abs(zc)), np.nan)


def rank_biserial(x, y):
    """Mann-Whitney rank-biserial correlation of x vs y (> 0: x larger) and two-sided p."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u / (len(x) * len(y)) - 1, p


def paired_rank_biserial(diff):
    """
    Matched-pairs rank-biserial correlation of paired differences and Wilcoxon signed-rank p, with zero differences
    kept (Pratt): all |differences| are ranked, zeros included, and r = (sum of positive ranks - sum of negative
    ranks) / sum of all ranks. Pairs without a difference (both 0, common for lowly expressed loci) thus pull r
    toward 0 instead of being ignored. NaN if every difference is 0.
    """
    d = np.asarray(diff, float)
    d = d[np.isfinite(d)]
    if not (d != 0).any():
        return np.nan, np.nan
    ranks = stats.rankdata(np.abs(d))
    r = (ranks[d > 0].sum() - ranks[d < 0].sum()) / ranks.sum()
    return r, stats.wilcoxon(d, zero_method="pratt").pvalue


def tcga_groups(samples):
    """Primary tumours and adjacent normals (one each per patient, with a `patient` barcode column) and GTEx liver."""
    t = samples[samples["source"] == "tcga"].copy()
    t["patient"] = t["tcga_barcode"].str[:12]
    one_per_patient = lambda kind: (t[t["cgc_sample_sample_type"] == kind].sort_values("tcga_barcode")
                                    .drop_duplicates("patient"))
    gtex_liver = samples[(samples["source"] == "gtex") & (samples["project"] == "LIVER")]
    return one_per_patient("Primary Tumor"), one_per_patient("Solid Tissue Normal"), gtex_liver


def d1_hcc(expr, samples, ncov=None, min_detected=MIN_DETECTED):
    """
    D1 HCC per feature (module docstring): paired tumour - normal (log2FC_paired, r_paired = D1_HCC, p, FDR) and,
    reported only, tumours vs GTEx liver (log2FC_vs_gtex, auc_vs_gtex, r_vs_gtex, p, FDR), with the group sizes.
    `ncov` (linear scale): features detected (ncov > 0) in < `min_detected` of the paired samples are not tested
    (paired statistics missing).
    """
    tum, nor, liver = tcga_groups(samples)
    pairs = tum.merge(nor, on="patient", suffixes=("_t", "_n"))
    e_tum = expr[tum["sample_id"]]
    e_gtex = expr[liver["sample_id"]]
    diff = expr[pairs["sample_id_t"]].to_numpy() - expr[pairs["sample_id_n"]].to_numpy()
    paired = list(pairs["sample_id_t"]) + list(pairs["sample_id_n"])
    detected = ((ncov[paired] > 0).mean(axis=1) >= min_detected).reindex(expr.index).to_numpy() \
        if ncov is not None else np.ones(len(expr), bool)
    rows = []
    for i, fid in enumerate(expr.index):
        r_p, p_p = paired_rank_biserial(diff[i]) if detected[i] else (np.nan, np.nan)
        r_g, p_g = rank_biserial(e_tum.iloc[i], e_gtex.iloc[i])
        rows.append((fid, diff[i].mean(), r_p, p_p,
                     e_tum.iloc[i].mean() - e_gtex.iloc[i].mean(), (r_g + 1) / 2, r_g, p_g))
    out = pd.DataFrame(rows, columns=["feature_id", "log2FC_paired", "r_paired", "p_paired", "log2FC_vs_gtex",
                                      "auc_vs_gtex", "r_vs_gtex", "p_vs_gtex"]).set_index("feature_id")
    out["fdr_paired"] = bh(out["p_paired"])
    out["fdr_vs_gtex"] = bh(out["p_vs_gtex"])
    out["n_pairs"] = len(pairs)
    out["n_tumour"] = len(tum)
    out["n_gtex_liver"] = len(liver)
    out["D1_HCC"] = out["r_paired"]
    return out


STAGE_LATE = ("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV", "Stage IVA", "Stage IVB")
STAGE_EARLY = ("Stage I", "Stage II")


def survival_table(samples, cdr, size_factors):
    """Primary tumours with OS/PFI endpoints, age, late-stage indicator and log size factor (patients without stage
    dropped)."""
    tum, _, _ = tcga_groups(samples)
    c = cdr[cdr["type"] == "LIHC"].copy()
    c = c[c["ajcc_pathologic_tumor_stage"].isin(STAGE_LATE + STAGE_EARLY)]
    c["late_stage"] = c["ajcc_pathologic_tumor_stage"].isin(STAGE_LATE).astype(float)
    c["age"] = pd.to_numeric(c["age_at_initial_pathologic_diagnosis"], errors="coerce")
    for col in ("OS", "OS.time", "PFI", "PFI.time"):
        c[col] = pd.to_numeric(c[col], errors="coerce")
    t = tum.merge(c, left_on="patient", right_on="bcr_patient_barcode")
    t["log_sf"] = np.log(t["sample_id"].map(size_factors))
    return t.dropna(subset=["age", "log_sf"])


COVARIATES = ["age", "late_stage", "log_sf"]


def _cox(df, endpoint):
    """Hazard ratio, Wald z and p of the standardized expression `x` (adjusted for COVARIATES) and the n used."""
    from lifelines import CoxPHFitter
    d = df[["x", *COVARIATES, endpoint, f"{endpoint}.time"]].dropna()
    d = d[d[f"{endpoint}.time"] > 0]
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(d, duration_col=f"{endpoint}.time", event_col=endpoint)
    s = cph.summary.loc["x"]
    return s["exp(coef)"], s["z"], s["p"], len(d)


def d2_survival(expr, samples, cdr, size_factors, min_expressed=MIN_DETECTED, ncov=None):
    """
    D2 per feature (module docstring): HR / z / p / n and BH-FDR per endpoint (OS, PFI) and D2 = mean z.
    `ncov` (linear scale) decides which features are expressed enough to test; `size_factors`: by sample_id. A model
    that does not converge leaves its columns missing (fit_error_<endpoint> says why).
    """
    from lifelines.exceptions import ConvergenceError
    tab = survival_table(samples, cdr, size_factors)
    e = expr[tab["sample_id"]]
    expressed = (ncov[tab["sample_id"]] > 0).mean(axis=1) if ncov is not None else pd.Series(1.0, index=e.index)
    rows = []
    for fid in e.index:
        rec = {"feature_id": fid, "frac_tumours_expressed": expressed[fid]}
        x = e.loc[fid].to_numpy()
        if expressed[fid] >= min_expressed and np.std(x) > 0:
            df = tab[[*COVARIATES, "OS", "OS.time", "PFI", "PFI.time"]].assign(x=(x - x.mean()) / x.std())
            for ep in ("OS", "PFI"):
                try:
                    rec[f"HR_{ep}"], rec[f"z_{ep}"], rec[f"p_{ep}"], rec[f"n_{ep}"] = _cox(df, ep)
                except (ConvergenceError, np.linalg.LinAlgError) as exc:
                    rec[f"fit_error_{ep}"] = type(exc).__name__
        rows.append(rec)
    out = pd.DataFrame(rows).set_index("feature_id")
    for ep in ("OS", "PFI"):
        for col in (f"HR_{ep}", f"z_{ep}", f"p_{ep}", f"n_{ep}"):
            if col not in out:
                out[col] = np.nan
        out[f"fdr_{ep}"] = bh(out[f"p_{ep}"])
    out["D2"] = out[["z_OS", "z_PFI"]].mean(axis=1)
    return out
