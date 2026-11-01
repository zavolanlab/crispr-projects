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
                           "hepg2_median_lfc": [-0.8, 0.1]})
    orfs = pd.DataFrame({"chrom": ["chr1", "chr1", "chr1"], "strand": ["+", "-", "+"],
                         "orf_blocks": ["203-230", "90-120", "100-229"]})
    m = screens.match(orfs, screen, "p")
    assert m["p_orfs"].tolist() == ["A", "B", ""]
    assert m["p_value"].tolist()[:2] == [2.0, 0.0] and np.isnan(m.loc[2, "p_value"])
    assert m.loc[0, "p_hepg2_lfc"] == -0.8
