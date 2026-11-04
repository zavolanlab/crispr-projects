import pandas as pd

from lincdisease import homology as h


def _hit(q, parent, pid, evalue=1e-20, qcov=90, bits=100, source="gc", s="M", **kw):
    d = dict(qseqid=q, sseqid=f"{source}|{pid}|{parent}", pident=90.0, length=50, qlen=60, slen=100, qstart=1,
             qend=50, sstart=1, send=50, evalue=evalue, bitscore=bits, qcovhsp=qcov, scovhsp=50.0,
             qseq_gapped=s, sseq_gapped=s)
    d.update(kw)
    return d


def test_homology_table_self_hits_and_thresholds():
    orfs = pd.DataFrame({"orf_id": ["o1", "o2", "o3"], "chrom": ["chr1"] * 3,
                         "orf_start0": [100, 5000, 9000], "orf_end1": [200, 5100, 9100]})
    spans = pd.DataFrame({"gene_name": ["SELF", "PAR", "WEAK"], "chrom": ["chr1", "chr2", "chr3"],
                          "start0": [150, 0, 0], "end1": [300, 10, 10]})
    hits = _as_hits(pd.DataFrame([
        _hit("o1", "SELF", "ENST1", bits=500),                                   # own locus: dropped
        _hit("o1", "PAR", "ENST2", bits=90),
        _hit("o2", "WEAK", "ENST3", evalue=1e-2),                                # E too high
        _hit("o3", "PAR", "P1", qcov=10, source="up"),                           # coverage too low
    ]))
    res = h.homology_table(orfs, hits, spans)
    assert res.loc["o1", "parent"] == "PAR" and res.loc["o1", "n_self_hits"] == 1
    assert res.loc[["o2", "o3"], "parent"].isna().all()


def _as_hits(df):
    sid = df["sseqid"].str.split("|", expand=True)
    return df.assign(source=sid[0].map({"gc": "gencode", "up": "uniprot"}), parent_id=sid[1], parent=sid[2])
