import numpy as np
import pandas as pd
import pytest

from lincdisease import qc


def _orfs():
    base = dict(chrom="chr1", strand="+", orf_antisense_to_CDS="", antisense_coding_exon_overlap="",
                orf_overlaps_small_RNA="", resolved_multimapper_dominated="False",
                possible_3UTR_extension_of="", possible_5prime_extension_of="", frac_psites_unique=0.9)
    rows = [
        dict(base, orf_id="clean", orf_start0=100, orf_end1=200),
        dict(base, orf_id="antisense", orf_start0=1000, orf_end1=1100, orf_antisense_to_CDS="GENE1"),
        dict(base, orf_id="paralog", orf_start0=5000, orf_end1=5100, frac_psites_unique=0.8),
        dict(base, orf_id="artefact", orf_start0=5200, orf_end1=5300, frac_psites_unique=0.05,
             resolved_multimapper_dominated="True"),
        dict(base, orf_id="artefact_with_peptide", orf_start0=5400, orf_end1=5500, frac_psites_unique=0.0),
    ]
    return pd.DataFrame(rows)


SEGDUPS = pd.DataFrame({"chrom": ["chr1", "chr1"], "start0": [4900, 5150], "end1": [6000, 5250],
                        "fracMatch": [0.95, 0.995]})


def test_interval_overlap_max_handles_nesting():
    o = _orfs()
    m = qc.interval_overlap_max(o, SEGDUPS, "fracMatch")
    assert np.isnan(m[0]) and np.isnan(m[1])
    assert m[2] == pytest.approx(0.95)
    assert m[3] == pytest.approx(0.995)          # nested region, larger value


def test_tiers():
    o = _orfs()
    res = qc.run_qc(o, segdups=SEGDUPS, copy_peptides=pd.Series({"artefact_with_peptide": True}))
    tier = dict(zip(o["orf_id"], res["confidence_tier"]))
    assert tier == {"clean": "clean", "antisense": "caution", "paralog": "expressed_paralog",
                    "artefact": "likely_artefact", "artefact_with_peptide": "expressed_paralog"}


def test_copy_attribution():
    o = _orfs()
    res = qc.run_qc(o, segdups=SEGDUPS)
    ca = dict(zip(o["orf_id"], res["copy_attribution"]))
    assert ca["clean"] == 1.0
    # paralog: expression 0.8, age (0.98 - 0.95) / 0.08 = 0.375
    assert ca["paralog"] == pytest.approx((0.8 + 0.375) / 2)
    # recent duplication with no unique reads: 0
    assert ca["artefact"] == pytest.approx((0.05 + 0.0) / 2)


def test_not_quantifiable_is_caution():
    o = _orfs().iloc[[0]]
    res = qc.run_qc(o, quantifiable=pd.Series({"clean": False}))
    assert res["confidence_tier"].iloc[0] == "caution"
