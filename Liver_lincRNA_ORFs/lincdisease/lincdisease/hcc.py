"""
Replication of D1 HCC (LOGIC.md, Step 3) in independent HCC cohorts with tumours and adjacent non-tumour liver of the
same patients, from recount3 SRA coverage quantified like TCGA/GTEx (same regions, `ncov`). Reported next to D1 HCC,
not scored.

Cohorts (patient and tissue from the recount3 SRA metadata):
  GSE77314  SRP068976: 50 HCC patients, tumour and adjacent normal (one run each)
  GSE77509  SRP069212: 20 HBV-HCC patients with portal vein tumour thrombus, primary tumour, adjacent normal and
            thrombus (several runs per sample; thrombi not used)
  GSE94660  SRP099053: 21 HBV-HCC patients, tumour and non-neoplastic liver
  GSE56545  SRP040998: HCC patients after liver transplantation, tumour and adjacent normal (single-end)

ncov is divided by median-of-ratios size factors of the recount3 gene sums, per cohort (expression.py: tumours lose
the hepatocyte transcripts that dominate liver RNA), and the runs of one sample are averaged. Per cohort and feature: matched-pairs rank-biserial correlation r of
log2(ncov + pseudocount), tumour minus normal of the same patient, and Wilcoxon signed-rank p (as D1 HCC paired);
features detected (ncov > 0) in < MIN_DETECTED of the cohort's samples are not tested there. Combined over cohorts:
hcc_rep_r = mean r weighted by the number of pairs, Stouffer z (signed z from the two-sided p, weights sqrt(pairs)),
p and BH-FDR over features. `hcc_replicated` (replicated): same sign as D1 HCC and combined FDR < 0.05; missing where
no cohort tested the feature.
"""
import re

import numpy as np
import pandas as pd
from scipy import stats

from . import expression, quant
from .chronic import parse_attributes
from .disease import bh, paired_rank_biserial

COHORTS = {"GSE77314": "SRP068976", "GSE77509": "SRP069212", "GSE94660": "SRP099053", "GSE56545": "SRP040998"}
MIN_DETECTED = 0.2


def sample_labels(cohort, title, attributes):
    """(patient, tissue 'tumour' / 'normal' / None) of one run from its experiment title and sample attributes."""
    if cohort == "GSE77314":                     # 'GSM2049057: Patient 1 Normal; Homo sapiens; RNA-Seq'
        m = re.search(r"Patient (\d+) (Normal|Tumor)", title)
        return (m.group(1), "tumour" if m.group(2) == "Tumor" else "normal") if m else (None, None)
    if cohort == "GSE77509":
        tissue = {"primary tumor": "tumour", "adjacent normal": "normal"}.get(attributes.get("tissue"))
        return attributes.get("patient id"), tissue
    if cohort == "GSE94660":                     # 'Tumor tissue of P105' / 'Non-neoplastic liver tissue of P105'
        m = re.search(r"(Tumor|Non-neoplastic liver) tissue of (P\d+)", title)
        return (m.group(2), "tumour" if m.group(1) == "Tumor" else "normal") if m else (None, None)
    if cohort == "GSE56545":     # source_name 'HCC1 Recurrent Tumor', "HCC-C11 Non-recurrent Tumor's Adjacent Normal"
        src = attributes.get("source_name", "")
        m = re.match(r"(HCC(?:-C)?\d+) ", src)
        return (m.group(1), "normal" if "Adjacent Normal" in src else "tumour") if m else (None, None)
    raise ValueError(cohort)


def manifest():
    """Runs of the replication cohorts (quant.manifest format) with cohort, sample (experiment), patient, tissue."""
    parts = []
    for cohort, project in COHORTS.items():
        m = quant.manifest("sra", project, ["sample_attributes", "experiment_title", "experiment_acc"])
        attrs = parse_attributes(m["sample_attributes"])
        lab = [sample_labels(cohort, str(t), a.dropna().to_dict())
               for t, (_, a) in zip(m["experiment_title"], attrs.iterrows())]
        m["cohort"] = cohort
        m["patient"] = [cohort + ":" + p if p else None for p, _ in lab]
        m["tissue"] = [t for _, t in lab]
        parts.append(m.dropna(subset=["patient", "tissue"]))
    man = pd.concat(parts, ignore_index=True)
    # keep patients with both tissues
    both = man.groupby("patient")["tissue"].nunique().eq(2)
    return man[man["patient"].map(both)].reset_index(drop=True)


def size_factors(genes, samples):
    """Size factors of the runs (expression.size_factors), per cohort; `genes`: protein-coding gene sums x runs."""
    s = samples.set_index("sample_id")
    return pd.concat([expression.size_factors(genes[genes.columns.intersection(g.index)])
                      for _, g in s.groupby("cohort")])


def replicated(rep, d1_hcc):
    """True / False per feature of `rep` (replication()): same sign as D1 HCC and combined FDR < 0.05; missing where
    no replication cohort tested the feature."""
    same = np.sign(rep["hcc_rep_r"]) == np.sign(d1_hcc.reindex(rep.index))
    return (same & (rep["hcc_rep_fdr"] < 0.05)).astype(object).where(rep["hcc_rep_fdr"].notna())


def replication(ncov, samples, pseudocount=0.05, min_detected=MIN_DETECTED):
    """
    ncov: features x runs (quant.normalize, scaled by size_factors); samples: manifest rows (sample_id, cohort, experiment_acc, patient,
    tissue). Returns per feature <cohort>_r / _p / _pairs and the combined hcc_rep_r, hcc_rep_z, hcc_rep_p,
    hcc_rep_fdr, hcc_rep_n_cohorts, hcc_rep_pairs.
    """
    s = samples.set_index("sample_id").loc[ncov.columns]
    per_sample = ncov.T.groupby(s["experiment_acc"]).mean().T          # runs of one sample averaged
    meta = s.drop_duplicates("experiment_acc").set_index("experiment_acc")
    e = np.log2(per_sample + pseudocount)
    out = pd.DataFrame(index=ncov.index)
    for cohort in COHORTS:
        m = meta[meta["cohort"] == cohort]
        if m.empty:
            continue
        tum = m[m["tissue"] == "tumour"].reset_index().drop_duplicates("patient").set_index("patient")
        nor = m[m["tissue"] == "normal"].reset_index().drop_duplicates("patient").set_index("patient")
        pts = tum.index.intersection(nor.index)
        diff = (e[tum.loc[pts, "experiment_acc"]].to_numpy() - e[nor.loc[pts, "experiment_acc"]].to_numpy())
        detected = (per_sample[m.index] > 0).mean(axis=1) >= min_detected
        res = [paired_rank_biserial(d) if ok else (np.nan, np.nan) for d, ok in zip(diff, detected)]
        out[f"{cohort}_r"] = [r for r, _ in res]
        out[f"{cohort}_p"] = [p for _, p in res]
        out[f"{cohort}_pairs"] = np.where(detected, len(pts), 0)
    cohorts = [c for c in COHORTS if f"{c}_r" in out]
    r = out[[f"{c}_r" for c in cohorts]].to_numpy(float)
    p = out[[f"{c}_p" for c in cohorts]].to_numpy(float)
    n = out[[f"{c}_pairs" for c in cohorts]].to_numpy(float)
    ok = np.isfinite(r) & (n > 0)
    w = np.where(ok, n, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["hcc_rep_r"] = np.where(ok.any(axis=1), (np.where(ok, r, 0) * w).sum(axis=1) / w.sum(axis=1), np.nan)
        z = np.sign(r) * stats.norm.isf(np.clip(p, 1e-300, 1) / 2)
        sw = np.sqrt(w)
        zc = (np.where(ok, z, 0) * sw).sum(axis=1) / np.sqrt((sw ** 2).sum(axis=1))
    out["hcc_rep_z"] = np.where(ok.any(axis=1), zc, np.nan)
    out["hcc_rep_p"] = np.where(np.isfinite(out["hcc_rep_z"]), 2 * stats.norm.sf(np.abs(out["hcc_rep_z"])), np.nan)
    out["hcc_rep_fdr"] = bh(out["hcc_rep_p"])
    out["hcc_rep_n_cohorts"] = ok.sum(axis=1)
    out["hcc_rep_pairs"] = w.sum(axis=1)
    return out
