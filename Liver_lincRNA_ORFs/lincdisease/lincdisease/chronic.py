"""
D1 chronic liver disease (LOGIC.md, Step 3): expression of loci and controls in recount3 SRA cohorts of
MASLD/MASH and viral hepatitis liver biopsies, quantified like TCGA/GTEx (same regions, `ncov`).

Cohorts (phenotypes from the recount3 SRA `sample_attributes`):
  hoang   SRP197353 (GSE130970, Hoang et al. 2019): 78 NAFLD biopsies, fibrosis stage 0-4 and NAFLD activity score
  suppli  SRP186450 (GSE126848, Suppli et al. 2019): 14 normal-weight, 12 obese, 15 NAFL, 16 NASH
  hcv     SRP078450 (GSE84346, Boldanova et al. 2017): untreated biopsies, 22 HCV-infected, 6 non-HCV

Per cohort, a signed statistic in [-1, 1] on log2(ncov + pseudocount) (> 0 = higher with disease):
  hoang   Spearman rho vs fibrosis stage (rho vs NAS reported)
  suppli  rank-biserial r of NAFL + NASH vs normal-weight + obese (Mann-Whitney; r = 2 AUC - 1)
  hcv     rank-biserial r of HCV vs non-HCV
D1_chronic = weighted mean of the available cohort statistics (fibrosis 0.5, suppli 0.25, hcv 0.25); the combined
p-value is Stouffer's (z from the two-sided p with the sign of the statistic, weights sqrt(n)), BH-FDR over features.
"""
import numpy as np
import pandas as pd
from scipy import stats

from . import quant
from .disease import bh

COHORTS = {"hoang": "SRP197353", "suppli": "SRP186450", "hcv": "SRP078450"}
WEIGHTS = {"hoang": 0.5, "suppli": 0.25, "hcv": 0.25}
MIN_DETECTED = 0.2      # features detected (ncov > 0) in fewer samples of a cohort are not tested


def parse_attributes(values):
    """recount3 `sample_attributes` ('key;;value|key;;value') -> DataFrame of attributes."""
    rows = []
    for v in values:
        d = {}
        for item in (v.split("|") if isinstance(v, str) else []):
            k, _, val = item.partition(";;")
            d[k.strip()] = val.strip()
        rows.append(d)
    return pd.DataFrame(rows)


def phenotypes(man):
    """Manifest rows (one recount3 SRA project each) -> cohort, fibrosis_stage, nas, group (case/control/NaN)."""
    a = parse_attributes(man["sample_attributes"]).set_index(man.index)
    out = pd.DataFrame(index=man.index)
    out["cohort"] = man["project"].map({v: k for k, v in COHORTS.items()})
    get = lambda col: a[col] if col in a else pd.Series(np.nan, index=a.index)
    out["fibrosis_stage"] = pd.to_numeric(get("fibrosis stage"), errors="coerce")
    out["nas"] = pd.to_numeric(get("nafld activity score"), errors="coerce")
    disease = get("disease").str.lower()
    infection = get("infection status")
    treated = get("treatment").fillna("None") != "None"
    group = pd.Series(np.nan, index=man.index, dtype=object)
    s = out["cohort"] == "suppli"
    group[s & disease.isin(["nafld", "nash"])] = "case"
    group[s & disease.isin(["healthy", "obese"])] = "control"
    h = (out["cohort"] == "hcv") & ~treated
    group[h & (infection == "HCV")] = "case"
    group[h & (infection == "non-HCV")] = "control"
    out["group"] = group
    return out


def manifest():
    """Samples of the chronic-disease cohorts (quant.manifest format) with their phenotypes."""
    parts = [quant.manifest("sra", p, ["sample_attributes", "sample_title"]) for p in COHORTS.values()]
    man = pd.concat(parts, ignore_index=True)
    man = pd.concat([man, phenotypes(man)], axis=1)
    keep = man["fibrosis_stage"].notna() | man["group"].notna()     # drops IFN-treated HCV biopsies
    return man[keep].reset_index(drop=True)


def _spearman(expr, y):
    out = {}
    for fid, x in expr.iterrows():
        ok = x.notna() & y.notna()
        if ok.sum() < 10 or x[ok].nunique() < 2:
            continue
        r, p = stats.spearmanr(x[ok], y[ok])
        out[fid] = (r, p, ok.sum())
    return pd.DataFrame(out, index=["stat", "p", "n"]).T


def _rank_biserial(expr, case, control):
    out = {}
    for fid, x in expr.iterrows():
        a, b = x[case].dropna(), x[control].dropna()
        if len(a) < 3 or len(b) < 3 or pd.concat([a, b]).nunique() < 2:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        out[fid] = (2 * u / (len(a) * len(b)) - 1, p, len(a) + len(b))
    return pd.DataFrame(out, index=["stat", "p", "n"]).T


def cohort_stats(expr, ncov, pheno):
    """Per feature and cohort: signed statistic, p, n (columns <cohort>_r/_p/_n) plus hoang rho vs NAS."""
    res = {}
    for cohort in COHORTS:
        idx = pheno.index[pheno["cohort"] == cohort]
        if not len(idx):
            continue
        detected = (ncov[idx] > 0).mean(axis=1) >= MIN_DETECTED
        e = expr.loc[detected, idx]
        if cohort == "hoang":
            t = _spearman(e, pheno.loc[idx, "fibrosis_stage"])
            nas = _spearman(e, pheno.loc[idx, "nas"])
            res["hoang_nas_rho"], res["hoang_nas_p"] = nas["stat"], nas["p"]
        else:
            g = pheno.loc[idx, "group"]
            t = _rank_biserial(e, g.index[g == "case"], g.index[g == "control"])
        res[f"{cohort}_r"], res[f"{cohort}_p"], res[f"{cohort}_n"] = t["stat"], t["p"], t["n"]
    return pd.DataFrame(res).reindex(expr.index)


def combine(cs, weights=WEIGHTS):
    """D1_chronic (weighted mean of cohort statistics), Stouffer p and BH-FDR."""
    cohorts = [c for c in weights if f"{c}_r" in cs]
    r = cs[[f"{c}_r" for c in cohorts]].to_numpy(float)
    p = cs[[f"{c}_p" for c in cohorts]].to_numpy(float)
    n = cs[[f"{c}_n" for c in cohorts]].to_numpy(float)
    w = np.array([weights[c] for c in cohorts])
    ok = np.isfinite(r)
    wsum = (ok * w).sum(axis=1)
    out = pd.DataFrame(index=cs.index)
    out["D1_chronic"] = np.where(wsum > 0, np.nansum(np.where(ok, r * w, 0), axis=1) / np.where(wsum > 0, wsum, 1), np.nan)
    z = np.sign(r) * stats.norm.isf(np.clip(p, 1e-300, 1) / 2)
    sw = np.where(ok, np.sqrt(n), 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        zc = np.nansum(np.where(ok, z * sw, 0), axis=1) / np.sqrt((sw ** 2).sum(axis=1))
    zc = np.where(ok.any(axis=1), zc, np.nan)
    out["chronic_z"] = zc
    out["chronic_p"] = np.where(np.isfinite(zc), 2 * stats.norm.sf(np.abs(zc)), np.nan)
    out["chronic_fdr"] = bh(out["chronic_p"])
    out["chronic_n_cohorts"] = ok.sum(axis=1)
    return out


def d1_chronic(ncov, pheno, pseudocount=0.05):
    """ncov: features x samples (cohort samples); pheno: phenotypes() indexed by sample_id."""
    ncov = ncov[pheno.index.intersection(ncov.columns)]
    pheno = pheno.loc[ncov.columns]
    expr = np.log2(ncov + pseudocount)
    cs = cohort_stats(expr, ncov, pheno)
    return pd.concat([combine(cs), cs], axis=1)
