"""
Optional fitted weights of the disease score (LOGIC.md, "Weight fitting").

Labels: established liver-disease lncRNAs = background loci with >= MIN_PAPERS non-retracted liver papers (D5
literature counts); negatives = background loci with an HGNC-style symbol and no liver paper (loci named only by
an Ensembl or StringTie id cannot have papers, so they are not negatives); loci in between are not labelled.
D5 itself is not a predictor, so the labels are not re-used as features.

Model: logistic regression of the label on the component percentiles (D1, D2, D3, D4, D6; missing = 0.5, centred
at 0.5) with non-negative coefficients and an L2 penalty, fitted by L-BFGS-B. The penalty is chosen by stratified
k-fold cross-validated AUC; the coefficients refitted on all labelled loci and scaled to sum 1 are the weights.
The cross-validated AUC is reported next to the AUC of the default weights on the same loci (no fitting, so not
optimistic).
"""
import numpy as np
import pandas as pd
from scipy import optimize
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from . import integrate, scoring

COMPONENTS = ("D1", "D2", "D3", "D4", "D6")
MIN_PAPERS = 5
LAMBDAS = (0.01, 0.1, 1.0, 10.0)


def labels(d5, is_background, min_papers=MIN_PAPERS):
    """1 / 0 / NaN per feature from the D5 table (symbol, n_liver, n_retracted)."""
    papers = d5["n_liver"] - d5["n_retracted"].fillna(0)
    named = ~d5["symbol"].astype(str).str.match(r"^(MSTRG\.|ENSG\d)")
    y = pd.Series(np.nan, index=d5.index)
    y[papers >= min_papers] = 1.0
    y[(d5["n_liver"] == 0) & named] = 0.0
    return y[is_background.reindex(y.index).fillna(False).astype(bool)]


def _design(comp, components):
    return comp.reindex(columns=list(components)).astype(float).fillna(scoring.NEUTRAL).to_numpy() - scoring.NEUTRAL


def fit_logistic_nonneg(X, y, lam):
    """Intercept and non-negative coefficients minimising the mean log-loss + lam * |w|^2."""
    def loss(theta):
        b, w = theta[0], theta[1:]
        z = b + X @ w
        nll = np.mean(np.logaddexp(0, z) - y * z)
        p = 1 / (1 + np.exp(-z))
        grad = np.r_[np.mean(p - y), X.T @ (p - y) / len(y) + 2 * lam * w]
        return nll + lam * w @ w, grad
    theta0 = np.r_[np.log(y.mean() / (1 - y.mean())), np.zeros(X.shape[1])]
    res = optimize.minimize(loss, theta0, jac=True, method="L-BFGS-B",
                            bounds=[(None, None)] + [(0, None)] * X.shape[1])
    return res.x[0], res.x[1:]


def fit(comp, y, components=COMPONENTS, lambdas=LAMBDAS, n_folds=5, seed=0):
    """Fitted weights (Series summing to 1, NaN if all zero) and a summary dict."""
    lab = y.dropna()
    X, t = _design(comp.loc[lab.index], components), lab.to_numpy()
    folds = list(StratifiedKFold(n_folds, shuffle=True, random_state=seed).split(X, t))
    cv = {}
    for lam in lambdas:
        pred = np.empty(len(t))
        for tr, te in folds:
            b, w = fit_logistic_nonneg(X[tr], t[tr], lam)
            pred[te] = X[te] @ w if w.any() else 0.0
        cv[lam] = roc_auc_score(t, pred) if np.ptp(pred) > 0 else 0.5
    best = max(cv, key=cv.get)
    _, w = fit_logistic_nonneg(X, t, best)
    weights = pd.Series(w / w.sum() if w.sum() > 0 else np.nan, index=list(components))
    default = {k: v for k, v in integrate.DISEASE_WEIGHTS.items() if k in components}
    default_score = scoring.composite(comp.loc[lab.index].reindex(columns=list(components)), default)["score"]
    summary = {"n_pos": int(t.sum()), "n_neg": int((1 - t).sum()), "lambda": best, "cv_auc": cv[best],
               "default_auc": roc_auc_score(t, default_score), **{f"cv_auc_lambda_{k}": v for k, v in cv.items()},
               # single-component AUCs (missing = 0.5) show which component carries the fit
               **{f"auc_{c}": roc_auc_score(t, X[:, i]) for i, c in enumerate(components)}}
    return weights, summary


def fitted_score(comp, weights):
    """Disease score with the fitted weights (same composite rule as the default score)."""
    w = weights.dropna()
    if w.empty or w.sum() <= 0:
        return pd.Series(np.nan, index=comp.index)
    return scoring.composite(comp.reindex(columns=list(w.index)), w.to_dict())["score"]
