import gzip

import numpy as np
import pandas as pd
import pytest

from lincdisease import homology as h


def test_translate_and_revcomp():
    assert h.translate("ATGGCCTAA") == "MA*"
    assert h.revcomp("ATGCN") == "NGCAT"


def test_syn_sites_known_codons():
    assert h.SYN_SITES["ATG"] == 0                       # Met: no synonymous change
    assert h.SYN_SITES["TGG"] == 0                       # Trp
    assert h.SYN_SITES["GCT"] == pytest.approx(1.0)      # 4-fold degenerate third position
    assert h.SYN_SITES["TTA"] == pytest.approx(2 / 3)    # Leu: TTG and CTA


def test_codon_diffs_paths():
    assert h._codon_diffs("GCT", "GCC") == (1.0, 0.0)
    assert h._codon_diffs("GCT", "ACT") == (0.0, 1.0)
    # CTT (Leu) -> TTG (Leu), two paths: via TTT (Phe: N,N) and via CTG (Leu: S,S)
    assert h._codon_diffs("CTT", "TTG") == (1.0, 1.0)


def test_nei_gojobori_synonymous_only():
    q = "GCT" * 30
    s = "GCC" * 15 + "GCT" * 15
    pairs = h.codon_pairs(q, s, "A" * 30, "A" * 30, 1, 1)
    r = h.nei_gojobori(pairs)
    assert r["n_codons"] == 30 and r["dN"] == 0 and r["dS"] > 0 and r["kaks"] == 0


def test_nei_gojobori_nan_without_ds_and_too_short():
    q = "GCT" * 30
    s = "ACT" * 5 + "GCT" * 25
    r = h.nei_gojobori(h.codon_pairs(q, s, "A" * 30, "A" * 30, 1, 1))
    assert r["dN"] > 0 and np.isnan(r["kaks"])            # dS = 0
    assert np.isnan(h.nei_gojobori(h.codon_pairs(q, s, "A" * 30, "A" * 30, 1, 1), min_codons=31)["kaks"])


def test_codon_pairs_follow_gaps():
    q = "ATGGCTAAA"          # M A K
    s = "ATGGGGGCTAAA"       # M G A K
    pairs = h.codon_pairs(q, s, "M-AK", "MGAK", 1, 1)
    assert pairs == [("ATG", "ATG"), ("GCT", "GCT"), ("AAA", "AAA")]


def test_orf_cds_blocks_and_strand(tmp_path):
    genome = tmp_path / "g.fa.gz"
    with gzip.open(genome, "wt") as fh:
        fh.write(">chr1 x\nAAATGGC\nCTAAGGG\n>chr2\nCCCTTAGGCCATCCC\n")
    orfs = pd.DataFrame({"orf_id": ["a", "b"], "chrom": ["chr1", "chr2"], "strand": ["+", "-"],
                         "orf_blocks": ["2-6;7-12", "3-12"]})
    cds = h.orf_cds(orfs, str(genome))
    assert cds["a"] == "ATGG" + "CTAAG"
    assert cds["b"] == h.revcomp("TTAGGCCAT")


def _hit(q, parent, pid, evalue=1e-20, qcov=90, bits=100, source="gc", s="M", **kw):
    d = dict(qseqid=q, sseqid=f"{source}|{pid}|{parent}", pident=90.0, length=50, qlen=60, slen=100, qstart=1,
             qend=50, sstart=1, send=50, evalue=evalue, bitscore=bits, qcovhsp=qcov, scovhsp=50.0,
             qseq_gapped=s, sseq_gapped=s)
    d.update(kw)
    return d


def test_homology_table_self_hits_thresholds_and_kaks():
    orfs = pd.DataFrame({"orf_id": ["o1", "o2", "o3"], "chrom": ["chr1"] * 3,
                         "orf_start0": [100, 5000, 9000], "orf_end1": [200, 5100, 9100]})
    spans = pd.DataFrame({"gene_name": ["SELF", "PAR", "WEAK"], "chrom": ["chr1", "chr2", "chr3"],
                          "start0": [150, 0, 0], "end1": [300, 10, 10]})
    prot = "M" + "A" * 29
    q_cds = "ATG" + "GCT" * 29 + "TAA"
    p_cds = "ATG" + "GCC" * 10 + "GCT" * 19 + "TAA"
    hits = _as_hits(pd.DataFrame([
        _hit("o1", "SELF", "ENST1", bits=500),                                   # own locus: dropped
        _hit("o1", "PAR", "ENST2", bits=90, s=prot, qend=30, send=30),
        _hit("o2", "WEAK", "ENST3", evalue=1e-2),                                # E too high
        _hit("o3", "PAR", "P1", qcov=10, source="up"),                           # coverage too low
    ]))
    res = h.homology_table(orfs, hits, spans, cds=pd.Series({"o1": q_cds}), parent_cds={"ENST2": p_cds})
    assert res.loc["o1", "parent"] == "PAR" and res.loc["o1", "n_self_hits"] == 1
    assert res.loc[["o2", "o3"], "parent"].isna().all()
    assert res.loc["o1", "kaks_n_codons"] == 30
    assert res.loc["o1", "kaks"] == 0 and res.loc["o1", "dS"] > 0


def _as_hits(df):
    sid = df["sseqid"].str.split("|", expand=True)
    return df.assign(source=sid[0].map({"gc": "gencode", "up": "uniprot"}), parent_id=sid[1], parent=sid[2])


def test_check_translation():
    orfs = pd.DataFrame({"orf_id": ["a", "b"], "aa_seq": ["MA", "MK"]})
    ok = h.check_translation(orfs, pd.Series({"a": "ATGGCTTAA", "b": "ATGGCTTAA"}))
    assert ok["a"] and not ok["b"]


def test_cds_matches_masked_residues():
    s = "ATG" + "AAA" * 5 + "GCT" + "TAA"
    assert h.cds_matches(s, "MKXX-XKA", 1)
    assert not h.cds_matches(s, "MKKKKKQ", 1)
    assert not h.cds_matches(s, "MKKKKKAAAA", 1)
