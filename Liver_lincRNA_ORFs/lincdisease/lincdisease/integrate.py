"""
Locus-level disease score (LOGIC.md, Step 3): component statistics -> percentiles against the background
loci -> weighted composite shrunk toward 0.5 for missing components.

Background = all lincRNA loci of the lincorfz run; control genes ('ctrl:<gene>') are scored against it but
are not part of it. Signed components (D1, D2, D3) are ranked on magnitude, so a locus strongly down in disease
scores as high as one strongly up; their signs give the direction (> 0: up in tumour / with disease, worse
survival, expression following a copy gain), summarized in disease_direction, hcc_direction and chronic_direction.
"""
import numpy as np
import pandas as pd

from . import scoring

DISEASE_WEIGHTS = {"D1": 0.25, "D2": 0.25, "D3": 0.20, "D4": 0.10, "D5": 0.10, "D6": 0.10}
D1_WEIGHTS = {"D1_HCC": 0.7, "D1_chronic": 0.3}
SIGNED = ("D1_HCC", "D1_chronic", "D2", "D3")
DIRECTION_MIN_PERCENTILE = 0.75
# the two scores reported next to disease_score (fixed in advance like the default weights); D6 (liver context)
# is shared, D4/D5 are split by trait class / paper topic
HCC_WEIGHTS = {"D1_HCC": 0.25, "D2": 0.25, "D3": 0.20, "D4_hcc": 0.10, "D5_hcc": 0.10, "D6": 0.10}
CHRONIC_WEIGHTS = {"D1_chronic": 0.50, "D4_chronic": 0.20, "D5_chronic": 0.15, "D6": 0.15}


def _gtex_medians(expr, samples, exclude=()):
    """features x GTEx tissues: median expression per tissue (recount3 project)."""
    gtex = samples[(samples["source"] == "gtex") & ~samples["project"].isin(exclude)]
    return pd.DataFrame({t: expr[g["sample_id"]].median(axis=1) for t, g in gtex.groupby("project")})


def liver_context(expr, samples, cell_lines):
    """
    D6 statistics per feature: median log expression in GTEx liver, liver specificity (liver median minus the
    largest median of the other GTEx tissues) and max(log2 TPM HepG2, Huh7) from lincorfz.
    """
    med = _gtex_medians(expr, samples)
    out = pd.DataFrame(index=expr.index)
    out["gtex_liver_median"] = med["LIVER"]
    out["liver_specificity"] = med["LIVER"] - med.drop(columns="LIVER").max(axis=1)
    out["cell_line_log2TPM"] = cell_lines.reindex(expr.index)
    return out


def tumour_specificity(expr, samples, exclude=("TESTIS",)):
    """Median in TCGA-LIHC primary tumours minus the largest median over GTEx tissues (testis excluded)."""
    tum = samples[(samples["source"] == "tcga") & (samples["cgc_sample_sample_type"] == "Primary Tumor")]
    return expr[tum["sample_id"]].median(axis=1) - _gtex_medians(expr, samples, exclude).max(axis=1)


def disease_score(stats, is_background, weights=DISEASE_WEIGHTS, d1_weights=D1_WEIGHTS):
    """
    `stats`: features x statistics (D1_HCC, D1_chronic, D2, D3, D4, D5, their hcc / chronic parts and the D6
    parts). Returns the percentiles (pct_), the components (comp_, with the sub-composites D1 and D6),
    disease_score / _coverage / _score_noD5, hcc_score / _coverage, chronic_score / _coverage and the directions.
    """
    bg = stats[is_background]
    cols = list(stats.columns)
    pct = scoring.component_percentiles(stats, bg, cols, signed=[c for c in SIGNED if c in cols])
    sign_cols = [c for c in pct.columns if c.endswith("_sign")]
    signs = pct[sign_cols].rename(columns=lambda c: c[:-5])
    pct = pct.drop(columns=sign_cols)

    comp = pd.DataFrame(index=stats.index)
    part_cov = {}     # coverage of the sub-composites, so that a missing part of D1 / D6 counts against coverage
    d1 = scoring.composite(pct.reindex(columns=list(d1_weights)), d1_weights)
    comp["D1"] = d1["score"].where(d1["coverage"] > 0)
    part_cov["D1"] = d1["coverage"]
    d6_parts = [c for c in ("gtex_liver_median", "liver_specificity", "cell_line_log2TPM") if c in pct]
    if d6_parts:
        d6 = scoring.composite(pct[d6_parts], dict.fromkeys(d6_parts, 1.0))
        comp["D6"] = d6["score"].where(d6["coverage"] > 0)
        part_cov["D6"] = d6["coverage"]
    # single-statistic components, and the parts that the HCC and chronic scores use directly
    for c in ("D2", "D3", "D4", "D5", "D1_HCC", "D1_chronic", "D4_hcc", "D4_chronic", "D5_hcc", "D5_chronic"):
        comp[c] = pct[c] if c in pct else np.nan

    have = comp.notna().astype(float)
    for c, cv in part_cov.items():
        have[c] = cv

    def coverage(w):
        """Fraction of the weight with data, where D1 and D6 count with the coverage of their sub-composites (a
        locus with D1 HCC but not D1 chronic has D1 coverage 0.7); scoring.composite's own coverage would count
        them as fully present."""
        w = pd.Series(w, dtype=float)
        return (have.reindex(columns=w.index).fillna(0) * w).sum(axis=1) / w.sum()

    res = scoring.composite(comp, weights)
    no5 = {k: v for k, v in weights.items() if k != "D5"}
    out = pd.concat([pct.add_prefix("pct_"), comp.add_prefix("comp_")], axis=1)
    out["disease_score"] = res["score"]
    out["disease_coverage"] = coverage(weights)
    out["disease_score_noD5"] = scoring.composite(comp, no5)["score"]
    for name, w in (("hcc", HCC_WEIGHTS), ("chronic", CHRONIC_WEIGHTS)):
        r = scoring.composite(comp, w)
        out[f"{name}_score"], out[f"{name}_coverage"] = r["score"], coverage(w)

    def direction(cols, labels=("oncogene-like", "suppressor-like")):
        cols = [c for c in cols if c in signs]
        return scoring.consensus_direction(pct[cols], signs[cols], DIRECTION_MIN_PERCENTILE, labels=labels)

    out["disease_direction"] = direction(SIGNED)
    out["hcc_direction"] = direction(("D1_HCC", "D2", "D3"))
    if "D1_chronic" in signs:
        out["chronic_direction"] = direction(["D1_chronic"], labels=("up", "down"))
    return out
