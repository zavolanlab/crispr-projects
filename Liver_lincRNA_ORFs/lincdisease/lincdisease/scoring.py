"""
Scoring arithmetic shared by all disease, actionability and mechanism scores (LOGIC.md, "General rules").

* every component is a statistic turned into an empirical percentile among the background loci/ORFs
  (mid-rank: ties count one half); signed statistics are ranked on their magnitude and the sign is kept
  separately as a direction; background members can carry weights (ORFs: 1 / number of ORFs of their locus, so
  that every locus counts once);
* a composite score is the weighted mean of its component percentiles where missing components contribute
  0.5 (no evidence either way), so the score is shrunk toward 0.5 in proportion to the missing weight;
  `coverage` is the fraction of the weight with data;
* a mechanism class score is the composite of its feature percentiles with equal weights; classes are compared
  on rank-based standardized scores (normal quantile of the class score's percentile among the background).
"""
import numpy as np
import pandas as pd
from scipy import stats

NEUTRAL = 0.5


def percentile_vs_background(values, background, weights=None):
    """Mid-rank empirical percentile of each value among the non-missing background values (NaN stays NaN).
    `weights` (aligned with `background`): weight of each background value (default 1)."""
    bg = pd.Series(background, dtype=float).reset_index(drop=True)
    w = (pd.Series(1.0, index=bg.index) if weights is None
         else pd.Series(np.asarray(weights, dtype=float), index=bg.index))
    keep = bg.notna() & w.gt(0)
    bg, w = bg[keep].to_numpy(), w[keep].to_numpy()
    x = np.asarray(pd.Series(values, dtype=float))
    out = np.full(x.shape, np.nan)
    if bg.size == 0:
        return out
    order = np.argsort(bg, kind="stable")
    bg, w = bg[order], w[order]
    cw = np.concatenate([[0.0], np.cumsum(w)])
    ok = ~np.isnan(x)
    below = cw[np.searchsorted(bg, x[ok], side="left")]
    upto = cw[np.searchsorted(bg, x[ok], side="right")]
    out[ok] = (below + 0.5 * (upto - below)) / cw[-1]
    return out


def signed_percentile_vs_background(values, background, weights=None):
    """Percentile of |value| among |background| and the sign (+1, -1, 0; NaN if missing)."""
    x = pd.Series(values, dtype=float)
    pct = percentile_vs_background(x.abs(), pd.Series(background, dtype=float).abs(), weights)
    return pct, np.sign(np.asarray(x))


def component_percentiles(table, background, columns, signed=(), weights=None):
    """
    Percentiles of `columns` of `table` against the same columns of `background` (rows weighted by `weights`).
    Columns listed in `signed` are ranked on magnitude; their signs are returned as '<col>_sign'.
    """
    out = pd.DataFrame(index=table.index)
    for col in columns:
        if col in signed:
            pct, sign = signed_percentile_vs_background(table[col], background[col], weights)
            out[col] = pct
            out[f"{col}_sign"] = sign
        else:
            out[col] = percentile_vs_background(table[col], background[col], weights)
    return out


def composite(percentiles, weights):
    """
    Weighted mean of component percentiles with missing components set to 0.5.
    Returns a DataFrame with 'score' and 'coverage' (fraction of the total weight with data).
    """
    w = pd.Series(weights, dtype=float)
    if (w < 0).any() or w.sum() <= 0:
        raise ValueError("weights must be non-negative and not all zero")
    p = percentiles.reindex(columns=w.index).astype(float)
    have = p.notna()
    score = (p.fillna(NEUTRAL) * w).sum(axis=1) / w.sum()
    coverage = (have * w).sum(axis=1) / w.sum()
    return pd.DataFrame({"score": score, "coverage": coverage}, index=percentiles.index)


def class_scores(features, classes):
    """
    Absolute evidence score per mechanism class (equal-weight composite of its feature percentiles)
    and its coverage. `classes` maps class name -> list of feature columns (features absent from the
    table count as missing).
    """
    scores, coverage = {}, {}
    for name, cols in classes.items():
        res = composite(features.reindex(columns=cols), dict.fromkeys(cols, 1.0))
        scores[name] = res["score"]
        coverage[name] = res["coverage"]
    return pd.DataFrame(scores), pd.DataFrame(coverage)


def locus_weights(locus_ids):
    """1 / number of rows of the same locus: every locus has total weight 1."""
    ids = pd.Series(np.asarray(locus_ids))
    return (1.0 / ids.map(ids.value_counts())).to_numpy()


def standardized(scores, background, coverage=None, weights=None):
    """
    Rank-based standardized class scores: the normal quantile of each class score's mid-rank percentile among the
    background (weighted by `weights`), per class. Unlike (score - mean) / SD, a class with few or binary features
    cannot reach extreme values by its spread alone: the top of every class maps to the same quantiles. Rows whose
    class `coverage` is 0 (no data for the class) get no standardized score for it, and the background of the class
    is the rows with data.
    """
    background = np.asarray(background, dtype=bool)
    w = np.ones(len(scores)) if weights is None else np.asarray(weights, dtype=float)
    out = pd.DataFrame(index=scores.index, columns=scores.columns, dtype=float)
    for c in scores.columns:
        s = scores[c].astype(float)
        if coverage is not None:
            s = s.where(coverage[c] > 0)
        ref = background & s.notna().to_numpy()
        n = ref.sum()
        if not n:
            continue
        pct = percentile_vs_background(s, s[ref], w[ref])
        out[c] = stats.norm.ppf(np.clip(pct, 0.5 / n, 1 - 0.5 / n))
    return out


def dominant_class(z, margin):
    """Top class if its standardized score is above the background median (> 0) and exceeds the second one by
    >= `margin` (standard normal units), else 'none'. Classes without a score (no data) do not compete."""
    def pick(row):
        r = row.dropna().sort_values(ascending=False)
        if r.empty or r.iloc[0] <= 0:
            return "none"
        if len(r) == 1 or r.iloc[0] - r.iloc[1] >= margin:
            return r.index[0]
        return "none"
    return z.apply(pick, axis=1)


def consensus_direction(percentiles, signs, min_percentile, labels=("up", "down")):
    """
    Direction from the signs of informative components (percentile >= `min_percentile`):
    all informative signs positive -> labels[0], all negative -> labels[1], both -> 'inconsistent',
    none informative -> 'none'.
    """
    informative = percentiles.ge(min_percentile) & signs.notna() & signs.ne(0)
    pos = (informative & signs.gt(0)).any(axis=1)
    neg = (informative & signs.lt(0)).any(axis=1)
    out = pd.Series("none", index=percentiles.index, dtype=object)
    out[pos & ~neg] = labels[0]
    out[neg & ~pos] = labels[1]
    out[pos & neg] = "inconsistent"
    return out
