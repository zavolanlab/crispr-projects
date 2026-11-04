import numpy as np
import pandas as pd

from lincdisease import screens


def test_screen_value():
    screened = pd.Series([True, True, True, True, True, False])
    phenotype = pd.Series([1, 1, 1, 1, 0, 1])
    secondary = pd.Series(["Specific to ORF", "Specific to ORF and transcript subregions", "Specific to transcript",
                           np.nan, np.nan, "Specific to ORF"])
    v = screens.screen_value(screened, phenotype, secondary)
    assert v.tolist()[:5] == [2.0, 2.0, 0.0, 1.0, 0.0]
    assert np.isnan(v.iloc[5])


def test_match():
    screen = pd.DataFrame({"chrom": ["chr1", "chr1"], "strand": ["+", "-"], "name": ["A", "B"],
                           "blocks": [[(100, 110), (200, 230)], [(90, 150)]], "value": [2.0, 0.0],
                           "lfc": [-0.8, 0.1]})
    orfs = pd.DataFrame({"chrom": ["chr1", "chr1", "chr1"], "strand": ["+", "-", "+"],
                         "orf_blocks": ["203-230", "90-120", "100-229"]})
    m = screens.match(orfs, screen, "p", {"lfc": "min"})
    assert m["p_orfs"].tolist() == ["A", "B", ""]
    assert m["p_value"].tolist()[:2] == [2.0, 0.0] and np.isnan(m.loc[2, "p_value"])
    assert m.loc[0, "p_lfc"] == -0.8


def test_chen_hits():
    rng = np.random.default_rng(0)
    neg = rng.normal(0, 0.01, 500)
    guides = pd.DataFrame({"orf": ["hit"] * 10 + ["cut_only"] * 10,
                           "gamma": np.r_[np.full(10, -0.1), rng.normal(-0.015, 0.01, 10)]})
    r = screens.chen_hits(guides, neg)
    assert r.loc["hit", "hit"] and not r.loc["cut_only", "hit"]    # the cutting effect alone is not a hit
    assert abs(r.loc["hit", "phenotype"] + 0.1) < 1e-12


def test_lift_points(tmp_path):
    # hg19 chr1:100-200 -> hg38 chr1:1100-1200 (+); hg19 chr2:0-50 -> hg38 chr3 reverse (size 1000)
    chain = tmp_path / "c.chain"
    chain.write_text("chain 100 chr1 1000 + 100 200 chr1 2000 + 1100 1200 1\n100\n\n"
                     "chain 50 chr2 1000 + 0 50 chr3 1000 - 0 50 2\n50\n")
    ch = screens.read_chain(str(chain))
    r = screens.lift_points(ch, ["chr1", "chr2", "chr1"], [150, 10, 500], ["+", "+", "-"])
    assert r.loc[0].tolist() == ["chr1", 1150, "+"]
    assert r.loc[1].tolist() == ["chr3", 989, "-"]
    assert pd.isna(r.loc[2, "chrom"])


def test_combined_value():
    pv = pd.Series([0.0, np.nan, 1.0, np.nan])
    resolved = pd.Series([True, False, False, False])
    chen = pd.Series([1.0, 1.0, 0.0, np.nan])
    v = screens.combined_value(pv, resolved, chen)
    assert v.tolist()[:3] == [0.0, 1.0, 1.0] and np.isnan(v.iloc[3])
