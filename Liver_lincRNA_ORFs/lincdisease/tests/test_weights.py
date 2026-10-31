import numpy as np
import pandas as pd

from lincdisease import weights


def test_labels():
    d5 = pd.DataFrame({"symbol": ["HULC", "LINC00001", "MSTRG.5", "LINC00002"], "n_liver": [40, 0, 0, 2],
                       "n_retracted": [1, 0, 0, 0]}, index=["a", "b", "c", "d"])
    y = weights.labels(d5, pd.Series(True, index=d5.index))
    assert y["a"] == 1 and y["b"] == 0 and np.isnan(y["c"]) and np.isnan(y["d"])


def test_fit_recovers_informative_component():
    rng = np.random.default_rng(0)
    n = 600
    comp = pd.DataFrame(rng.uniform(0, 1, (n, 5)), columns=list(weights.COMPONENTS))
    p = 1 / (1 + np.exp(-(-2 + 6 * (comp["D2"] - 0.5))))
    y = pd.Series((rng.uniform(0, 1, n) < p).astype(float), index=comp.index)
    w, s = weights.fit(comp, y)
    assert w.idxmax() == "D2" and w["D2"] > 0.6 and abs(w.sum() - 1) < 1e-9
    assert s["cv_auc"] > 0.7 and (w >= 0).all()
    assert weights.fitted_score(comp, w).between(0, 1).all()
