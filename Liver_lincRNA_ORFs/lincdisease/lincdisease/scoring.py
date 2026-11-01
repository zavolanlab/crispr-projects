"""
Scoring arithmetic shared by all disease, actionability and mechanism scores (LOGIC.md, "General rules").

* every component is a statistic turned into an empirical percentile among the background loci/ORFs
  (mid-rank: ties count one half); signed statistics are ranked on their magnitude and the sign is kept
  separately as a direction;
* a composite score is the weighted mean of its component percentiles where missing components contribute
  0.5 (no evidence either way), so the score is shrunk toward 0.5 in proportion to the missing weight;
  `coverage` is the fraction of the weight with data;
* a mechanism class score is the composite of its feature percentiles with equal weights; its relative share
  is the score divided by the sum over classes.
"""
import numpy as np
import pandas as pd

NEUTRAL = 0.5


def percentile_vs_background(values, background):
    """Mid-rank empirical percentile of each value among the non-missing background values (NaN stays NaN)."""
    bg = np.sort(np.asarray(pd.Series(background, dtype=float).dropna()))
    x = np.asarray(pd.Series(values, dtype=float))
    out = np.full(x.shape, np.nan)
    if bg.size == 0:
        return out
    ok = ~np.isnan(x)
    below = np.searchsorted(bg, x[ok], side="left")
    upto = np.searchsorted(bg, x[ok], side="right")
    out[ok] = (below + 0.5 * (upto - below)) / bg.size
    return out


def signed_percentile_vs_background(values, background):
    """Percentile of |value| among |background| and the sign (+1, -1, 0; NaN if missing)."""
    x = pd.Series(values, dtype=float)
    pct = percentile_vs_background(x.abs(), pd.Series(background, dtype=float).abs())
    return pct, np.sign(np.asarray(x))


def component_percentiles(table, background, columns, signed=()):
    """
    Percentiles of `columns` of `table` against the same columns of `background`.
    Columns listed in `signed` are ranked on magnitude; their signs are returned as '<col>_sign'.
    """
    out = pd.DataFrame(index=table.index)
    for col in columns:
        if col in signed:
            pct, sign = signed_percentile_vs_background(table[col], background[col])
            out[col] = pct
            out[f"{col}_sign"] = sign
        else:
            out[col] = percentile_vs_background(table[col], background[col])
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


def standardized(scores, background):
    """Class scores in standard deviations from the background mean, per class, so that classes with few or
    coarse features (wider spread) do not win the comparison between classes by their spread alone."""
    ref = scores[background]
    return (scores - ref.mean()) / ref.std()


def dominant_class(z, margin):
    """Top class if its standardized score is above the background mean (> 0) and exceeds the second one by
    >= `margin` SD, else 'none'."""
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
