"""
Locus-level disease score (LOGIC.md, Step 3): component statistics -> percentiles against the background
loci -> weighted composite shrunk toward 0.5 for missing components.

Background = all lincRNA loci of the lincorfz run; control genes ('ctrl:<gene>') are scored against it but
are not part of it. Signed components (D1, D2, D3) are ranked on magnitude; their signs give the direction.
"""
import numpy as np
import pandas as pd

from . import scoring

DISEASE_WEIGHTS = {"D1": 0.25, "D2": 0.25, "D3": 0.20, "D4": 0.10, "D5": 0.10, "D6": 0.10}
D1_WEIGHTS = {"D1_HCC": 0.7, "D1_chronic": 0.3}
SIGNED = ("D1_HCC", "D1_chronic", "D2", "D3")
DIRECTION_MIN_PERCENTILE = 0.75


def liver_context(expr, samples, cell_lines):
    """
    D6 statistics per feature: median log expression in GTEx liver, liver specificity (liver median minus the
    largest median of the other GTEx tissues) and max(log2 TPM HepG2, Huh7) from lincorfz.
    """
    gtex = samples[samples["source"] == "gtex"]
    med = pd.DataFrame({t: expr[g["sample_id"]].median(axis=1) for t, g in gtex.groupby("project")})
    out = pd.DataFrame(index=expr.index)
    out["gtex_liver_median"] = med["LIVER"]
    out["liver_specificity"] = med["LIVER"] - med.drop(columns="LIVER").max(axis=1)
    out["cell_line_log2TPM"] = cell_lines.reindex(expr.index)
    return out


def tumour_specificity(expr, samples, exclude=("TESTIS",)):
    """Median in TCGA-LIHC primary tumours minus the largest median over GTEx tissues (testis excluded)."""
    tum = samples[(samples["source"] == "tcga") & (samples["cgc_sample_sample_type"] == "Primary Tumor")]
    gtex = samples[(samples["source"] == "gtex") & ~samples["project"].isin(exclude)]
    med = pd.DataFrame({t: expr[g["sample_id"]].median(axis=1) for t, g in gtex.groupby("project")})
    return expr[tum["sample_id"]].median(axis=1) - med.max(axis=1)


def disease_score(stats, is_background, weights=DISEASE_WEIGHTS, d1_weights=D1_WEIGHTS):
    """
    `stats`: features x statistics (D1_HCC, D1_chronic, D2, D3, D4, D5, and the D6 parts). Returns percentiles,
    signs, the sub-composites D1 and D6, disease_score, disease_coverage, disease_score_noD5 and direction.
    """
    bg = stats[is_background]
    cols = [c for c in stats.columns]
    pct = scoring.component_percentiles(stats, bg, cols, signed=[c for c in SIGNED if c in cols])
    sign_cols = [c for c in pct.columns if c.endswith("_sign")]
    signs = pct[sign_cols].rename(columns=lambda c: c[:-5])
    pct = pct.drop(columns=sign_cols)

    comp = pd.DataFrame(index=stats.index)
    d1 = scoring.composite(pct.reindex(columns=list(d1_weights)), d1_weights)
    comp["D1"] = d1["score"].where(d1["coverage"] > 0)
    d6_parts = [c for c in ("gtex_liver_median", "liver_specificity", "cell_line_log2TPM") if c in pct]
    if d6_parts:
        d6 = scoring.composite(pct[d6_parts], dict.fromkeys(d6_parts, 1.0))
        comp["D6"] = d6["score"].where(d6["coverage"] > 0)
    for c in ("D2", "D3", "D4", "D5"):
        comp[c] = pct[c] if c in pct else np.nan

    res = scoring.composite(comp, weights)
    no5 = {k: v for k, v in weights.items() if k != "D5"}
    out = pd.concat([pct.add_prefix("pct_"), comp.add_prefix("comp_")], axis=1)
    out["disease_score"] = res["score"]
    out["disease_coverage"] = res["coverage"]
    out["disease_score_noD5"] = scoring.composite(comp, no5)["score"]
    dir_cols = [c for c in SIGNED if c in signs]
    out["disease_direction"] = scoring.consensus_direction(pct[dir_cols], signs[dir_cols], DIRECTION_MIN_PERCENTILE,
                                                           labels=("oncogene-like", "suppressor-like"))
    # D3 sign is set so that > 0 means oncogene-like (CN gain / loss of methylation silencing raise expression)
    return out
