import numpy as np
import pandas as pd

from lincdisease import somatic


def test_tss_windows_strand():
    t = pd.DataFrame({"feature_id": ["a", "b"], "chrom": "chr1", "start0": [1000, 5000], "end1": [2000, 6000],
                      "strand": ["+", "-"]})
    w = somatic.tss_windows(t, upstream=100, downstream=10)
    assert list(w["start0"]) == [900, 5989] and list(w["end1"]) == [1011, 6100]


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
    cn_gain = rng.uniform(0, 1, n)            # mostly gained, dosage-driven
    beta = rng.uniform(0, 1, n)               # silencing by methylation
    expr = pd.DataFrame([cn_gain + rng.normal(0, 0.05, n), 1 - beta + rng.normal(0, 0.05, n)],
                        index=["onc", "sup"], columns=samples["sample_id"])
    cn = pd.DataFrame([cn_gain, np.zeros(n)], index=["onc", "sup"], columns=ids)
    normals = [f"TCGA-AA-{i:04d}-11A" for i in range(10)]
    b = pd.DataFrame(np.c_[np.vstack([np.full(n, 0.5), beta]), np.vstack([np.full(10, 0.5), np.full(10, 0.1)])],
                     index=["onc", "sup"], columns=ids + normals)
    res = somatic.d3(expr, samples, cn, b)
    assert res.loc["onc", "D3"] > 0.3 and res.loc["onc", "D3_source"] == "cn"
    assert res.loc["sup", "D3"] < 0 and res.loc["sup", "D3_source"] == "meth"


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


def test_dominant_isoforms():
    t = pd.DataFrame({"feature_id": ["a", "a", "ctrl:B", "ctrl:B"], "transcript_id": ["a1", "a2", "b1", "b2"],
                      "TPM_max": [1.0, 5.0, np.nan, np.nan], "canonical": [False, False, False, True]})
    assert list(somatic.dominant_isoforms(t)["transcript_id"]) == ["a2", "b2"]
