import json

import numpy as np
import pandas as pd

from lincdisease import table_html


def _tables():
    orfs = pd.DataFrame({"orf_id": ["a", "b", "c"], "locus_id": ["L1", "L1", "L2"], "locus_name": ["X", "X", "Y"],
                         "is_candidate": [True, False, False], "disease_score_noD5": [0.6, 0.9, np.nan],
                         "hcc_score": [0.6, 0.9, np.nan],
                         "chronic_score": [0.5, 0.4, 0.7], "pct_m1_orf_share": [0.1, 0.2, 0.3]})
    support = pd.DataFrame({"orf_id": ["a", "b", "c"], "crispr_screen_value": [np.nan, 2.0, 0.0]})
    mech = pd.DataFrame({"orf_id": ["a", "b", "c"], "pct_m1_orf_share": [0.1, 0.2, 0.3], "pct_m2_nmd": [1, 0, 0.5]})
    disease = pd.DataFrame({"stat_D3": [0.2, -0.1], "cn_rho": [0.3, 0.0]}, index=["L1", "L2"])
    peptides = pd.DataFrame({"n_specific_peptides": [1, 0, 0], "n_hla_specific": [0, 0, 2]}, index=["a", "b", "c"])
    hla = pd.DataFrame({"hla_any_n_binders": [3, 0, 1]}, index=["a", "b", "c"])
    return orfs, support, mech, disease, peptides, hla


def test_rows_union_and_sets():
    t = table_html.rows(*_tables())
    assert list(t["orf_id"]) == ["b", "a"]           # candidates + screen value >= 1, by disease score w/o D5
    assert t.set_index("orf_id")["set"].to_dict() == {"a": "candidate", "b": "screen"}
    assert t.set_index("orf_id").loc["b", "stat_D3"] == 0.2 and "pct_m2_nmd" in t


def test_build_orders_summary_before_detail_and_encodes_missing():
    t = table_html.rows(*_tables())
    page = table_html.build(t)
    meta = json.loads(page.split("const META = ")[1].split(";\n")[0])
    data = json.loads(page.split("const DATA = ")[1].split(";\n")[0])
    cols = [m["col"] for m in meta]
    assert cols.index("hcc_score") < cols.index("chronic_score") < cols.index("stat_D3") < cols.index("pct_m2_nmd")
    assert len(data) == 2 and None in data[1]          # NaN becomes null, not NaN (valid JSON)
