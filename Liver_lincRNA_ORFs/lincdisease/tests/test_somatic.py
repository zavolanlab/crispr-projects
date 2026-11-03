import numpy as np
import pandas as pd

from lincdisease import somatic


def test_cn_at_midpoint():
    spans = pd.DataFrame({"chrom": ["chr1", "chr1", "chr2"], "start0": [0, 300, 0], "end1": [100, 500, 10]},
                         index=["a", "b", "c"])
    seg = pd.DataFrame({"sample": "S-01A", "Chrom": ["1", "1"], "Start": [1, 201], "End": [200, 1000],
                        "value": [0.5, -0.4]})
    cn = somatic.cn_at(spans, seg)
    assert cn.loc["a", "S-01A"] == 0.5 and cn.loc["b", "S-01A"] == -0.4 and np.isnan(cn.loc["c", "S-01A"])


def test_d3_sign():
    rng = np.random.default_rng(0)
    n = 60
    ids = [f"TCGA-AA-{i:04d}-01A" for i in range(n)]
    samples = pd.DataFrame({"source": "tcga", "sample_id": [f"tcga:{i}" for i in range(n)],
                            "tcga_barcode": [b + "-11R" for b in ids]})
    gain = rng.uniform(0, 1, n)               # mostly gained, dosage-driven
    loss = -rng.uniform(0, 1, n)              # mostly lost, dosage-driven
    flat = rng.uniform(0, 1, n)               # gained, expression unrelated
    expr = pd.DataFrame([gain + rng.normal(0, 0.05, n), loss + rng.normal(0, 0.05, n), rng.normal(0, 1, n)],
                        index=["onc", "sup", "none"], columns=samples["sample_id"])
    cn = pd.DataFrame([gain, loss, flat], index=["onc", "sup", "none"], columns=ids)
    res = somatic.d3(expr, samples, cn)
    assert res.loc["onc", "D3"] > 0.3
    assert res.loc["sup", "D3"] < -0.3
    assert abs(res.loc["none", "D3"]) < 0.2


def test_partial_spearman_removes_confounder():
    rng = np.random.default_rng(1)
    z = rng.uniform(0, 1, 300)                           # purity drives both
    x, y = z + rng.normal(0, 0.1, 300), z + rng.normal(0, 0.1, 300)
    assert stats_rho(x, y) > 0.8
    r, p = somatic.partial_spearman(x, y, z)
    assert abs(r) < 0.15 and p > 0.01


def stats_rho(x, y):
    from scipy import stats
    return stats.spearmanr(x, y)[0]
