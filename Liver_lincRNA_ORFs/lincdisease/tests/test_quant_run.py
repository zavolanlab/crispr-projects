"""quant.run: resuming a sums file, column checks against the regions of the existing sums; gene sums checks."""
import sys
import types

import pandas as pd
import pytest

from lincdisease import __main__ as cli
from lincdisease import quant


class FakeBW:
    def chroms(self):
        return {"c1": 1000}

    def stats(self, c, s, e, type, exact):
        return [float(e - s)]          # coverage 1 everywhere

    def close(self):
        pass


@pytest.fixture
def fake_bigwig(monkeypatch):
    monkeypatch.setitem(sys.modules, "pyBigWig", types.SimpleNamespace(open=lambda path: FakeBW()))
    monkeypatch.setattr(quant, "_download", lambda url, dest: open(dest, "w").close())


def _regions(features=("a", "b")):
    rows = {"a": ("c1", 0, 50), "b": ("c1", 100, 110), "c": ("c1", 200, 220)}
    return pd.DataFrame([(f, *rows[f]) for f in features], columns=["feature_id", "chrom", "start0", "end1"])


def _samples(ids):
    return pd.DataFrame({"sample_id": ids, "url": "x"})


def _sums(d):
    return pd.read_csv(d / "sums.0.tsv", sep="\t")


def test_run_resumes_and_drops_incomplete_rows(tmp_path, fake_bigwig):
    quant.run(_samples(["s1", "s2"]), _regions(), tmp_path)
    with open(tmp_path / "sums.0.tsv", "a") as fh:
        fh.write("s3\t50")                         # cut short by a kill: no line end
    quant.run(_samples(["s1", "s2", "s3"]), _regions(), tmp_path)
    s = _sums(tmp_path)
    assert list(s["sample_id"]) == ["s1", "s2", "s3"]
    assert s.set_index("sample_id").loc["s3"].to_dict() == {"a": 50.0, "b": 10.0}


def test_run_keeps_the_columns_of_existing_sums(tmp_path, fake_bigwig):
    quant.run(_samples(["s1"]), _regions(("a", "b")), tmp_path)
    quant.run(_samples(["s1", "s2"]), _regions(("b",)), tmp_path)     # feature a no longer quantified
    s = _sums(tmp_path).set_index("sample_id")
    assert list(s.columns) == ["a", "b"] and s.loc["s2", "b"] == 10.0 and pd.isna(s.loc["s2", "a"])


def test_run_refuses_new_features_and_changed_regions(tmp_path, fake_bigwig):
    quant.run(_samples(["s1"]), _regions(("a", "b")), tmp_path)
    with pytest.raises(ValueError, match="not in the existing sums"):
        quant.run(_samples(["s1", "s2"]), _regions(("a", "b", "c")), tmp_path)
    moved = _regions(("a", "b"))
    moved.loc[1, "end1"] = 120
    with pytest.raises(ValueError, match="differ"):
        quant.run(_samples(["s1", "s2"]), moved, tmp_path)


def test_run_refuses_a_sums_file_with_other_columns(tmp_path, fake_bigwig):
    quant.run(_samples(["s1"]), _regions(("a", "b")), tmp_path)
    (tmp_path / "sums.0.tsv").write_text("sample_id\tb\ta\ns1\t10\t50\n")
    with pytest.raises(ValueError, match="columns differ"):
        quant.run(_samples(["s1", "s2"]), _regions(("a", "b")), tmp_path)


def test_gene_sums_raises_on_missing_samples(tmp_path, monkeypatch):
    path = tmp_path / "genes.tsv"
    path.write_text("##recount3\ngene_id\tr1\tr2\nENSG1.1\t10\t20\n")
    monkeypatch.setattr(cli, "fetch_gene_sums", lambda data_dir, source, project: str(path))
    samples = pd.DataFrame({"sample_id": ["sra:r1", "sra:r2", "sra:r3", "tcga:t1"],
                            "source": ["sra", "sra", "sra", "tcga"], "project": ["P1", "P1", "P1", "LIHC"],
                            "auc": 1e9})
    with pytest.raises(ValueError, match="1 samples"):
        cli.gene_sums(tmp_path, [("sra", "P1")], samples)
    genes = cli.gene_sums(tmp_path, [("sra", "P1")], samples[samples["sample_id"] != "sra:r3"])
    assert list(genes.columns) == ["sra:r1", "sra:r2"]           # tcga:t1 belongs to another project
