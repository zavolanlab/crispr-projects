import numpy as np
import pandas as pd
import pytest

from lincdisease import hla, peptide


def test_parse_blocks_and_stop_end():
    b = peptide.parse_blocks("100-110;200-230")
    assert b == [(100, 110), (200, 230)]
    assert peptide.stop_end(b, "+") == 230
    assert peptide.stop_end(b, "-") == 100


def test_same_stop(tmp_path):
    bed = tmp_path / "c.bed"
    bed.write_text("chr1\t100\t230\tc1orf1:G:lncRNA\t0\t+\t100\t230\t0\t2\t10,30\t0,100\n"
                   "chr1\t90\t150\tc1orf2:G:lncRNA\t0\t-\t90\t150\t0\t1\t60\t0\n")
    cat = peptide.read_bed12(str(bed))
    assert cat.loc[0, "blocks"] == [(100, 110), (200, 230)]
    orfs = pd.DataFrame({"chrom": ["chr1", "chr1", "chr1", "chr2"], "strand": ["+", "-", "+", "+"],
                         "orf_blocks": ["203-230", "90-120", "100-229", "100-230"]})
    # same stop with a different start; minus-strand stop at the lowest coordinate; stop shifted; other chrom
    assert peptide.same_stop(orfs, cat).tolist() == ["c1orf1", "c1orf2", "", ""]


def test_tm_segments():
    hydrophobic = "L" * 22
    polar = "DEKRSNQ" * 4
    assert peptide.tm_segments(polar + hydrophobic + polar) == 1
    assert peptide.tm_segments(polar + hydrophobic + polar + hydrophobic + polar) == 2
    assert peptide.tm_segments(polar * 3) == 0
    assert peptide.tm_segments("LLL") == 0


def test_signal_like():
    albumin = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQ"
    assert peptide.signal_like(albumin)
    assert not peptide.signal_like("MSTVHEILCKLSLEGDHSTPPSAYGSVKAYTNFDAERDALNIETAIKTKGVDEV")
    # hydrophobic stretch but no small residues at -3/-1 downstream
    assert not peptide.signal_like("M" + "KR" + "L" * 10 + "QQQQQQQQQQQQQQQQQ")


def test_sequence_features():
    orfs = pd.DataFrame({"aa_seq": ["MKWVTFISLLFLFSSAYSRGVFRRDAH*", "MDEK*"]})
    f = peptide.sequence_features(orfs)
    assert f["pep_length"].tolist() == [27, 4]
    assert f["pep_signal_like"].tolist() == [True, False]
    assert np.isnan(f.loc[1, "pep_max_hydropathy19"])


def test_read_domtblout(tmp_path):
    line = "orf1 - 120 7tm_1 PF00001.24 260 1e-30 100.0 0.1 1 1 1e-31 1e-30 99.0 0.1 10 250 5 110 4 111 0.95 -\n"
    p = tmp_path / "d.tsv"
    p.write_text("# header\n" + line)
    d = peptide.read_domtblout(str(p))
    assert d.loc[0, "key"] == "orf1" and d.loc[0, "pfam_name"] == "7tm_1" and d.loc[0, "pfam_acc"] == "PF00001"
    assert (d.loc[0, "env_from"], d.loc[0, "env_to"]) == (4, 111)


def test_hla_summarize():
    pred = pd.DataFrame({
        "sequence_name": ["a", "a", "a", "a", "b"],
        "sample_name": ["HepG2", "HepG2", "Huh7", "common:HLA-A*01:01", "HepG2"],
        "peptide": ["PEPTIDEA", "PEPTIDEB", "PEPTIDEA", "PEPTIDEC", "PEPTIDED"],
        "presentation_percentile": [0.5, 3.0, 1.0, 1.9, 10.0],
    })
    _, set_of = hla.samples()
    s = hla.summarize(pred, set_of)
    assert s.loc["a", "hla_HepG2_n_binders"] == 1
    assert s.loc["a", "hla_Huh7_n_binders"] == 1
    assert s.loc["a", "hla_common_n_binders"] == 1
    assert s.loc["a", "hla_any_n_binders"] == 2          # PEPTIDEA counted once
    assert s.loc["b", "hla_HepG2_n_binders"] == 0
    assert s.loc["b", "hla_HepG2_best_pct"] == pytest.approx(10.0)
    f = hla.finalize(s, ["a", "b", "short"])
    assert f.loc["short", "hla_any_n_binders"] == 0 and np.isnan(f.loc["short", "hla_HepG2_best_pct"])
