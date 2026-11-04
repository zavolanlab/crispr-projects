import numpy as np
import pandas as pd

from lincdisease import controls, report


def _disease(**ctrl_values):
    """disease.loci.tsv-like table: 10 background loci and the HULC control with the given values."""
    bg = pd.DataFrame({"disease_score": np.linspace(0, 1, 10), "is_control": False},
                      index=[f"L{i}" for i in range(10)])
    ctrl = pd.DataFrame({"is_control": True, **{k: [v] for k, v in ctrl_values.items()}}, index=["ctrl:HULC"])
    return pd.concat([bg, ctrl])


def _rows(table, score, gene="HULC"):
    return table[(table["score"] == score) & (table["control"] == gene)]


def test_control_report_pass_fail_and_not_evaluated():
    d = _disease(log2FC_paired=0.5, fdr_paired=0.01, disease_score=0.5, pct_D5=np.nan)
    mech = pd.DataFrame({**{f"M{i}_z": [0.0] for i in range(1, 7)}, **{f"M{i}_score": [0.5] for i in range(1, 7)}},
                        index=["GAS5"])
    t = report.control_report(d, mech)
    assert _rows(t, "D1_HCC")["pass"].item() is True
    assert _rows(t, "disease_score")["pass"].item() is False          # 0.5 < background 90th percentile
    assert pd.isna(_rows(t, "D5")["pass"].item())                     # value missing: not evaluated
    assert pd.isna(_rows(t, "D2", "BIRC5")["pass"].item())            # control not quantified
    assert pd.isna(_rows(t, "M1", "MTLN")["pass"].item())             # not a locus of the run
    assert "log2FC_paired > 0" in _rows(t, "D1_HCC")["expectation"].item()
    assert len(t) == len(controls.PANEL)


def test_control_genes_cover_panel_and_neighbours():
    genes = controls.control_genes()
    assert {"HULC", "GPC3", "HNF1A", "FOXA2"} <= set(genes) and genes == sorted(set(genes))
