import json
import re

import pandas as pd

from lincdisease import literature


def test_usable_alias():
    approved = {"TMEM99X"}
    assert not literature.usable_alias("ALIEN", approved)      # a word
    assert not literature.usable_alias("LEU2", approved)       # short, yeast gene
    assert not literature.usable_alias("TMEM99X", approved)    # another gene's symbol
    assert literature.usable_alias("LINC00657", approved)


def test_queries_are_exact_phrases():
    q = literature.build_queries("LINC01090")
    assert '(TITLE:"LINC01090" OR ABSTRACT:"LINC01090")' in q["n_liver"]
    assert "Retracted Publication" in q["n_retracted"]


def test_usable_alias_drops_generic_tokens():
    assert not literature.usable_alias("BETA-3", set())
    assert not literature.usable_alias("EEF-1", set())
    assert literature.usable_alias("MALAT-1", set())


def test_hgnc_aliases_drop_shared_and_stem_names(tmp_path):
    path = tmp_path / "hgnc.txt"
    pd.DataFrame({"symbol": ["EEF1B2", "EEF1B2P2", "GAS5", "SNHG20", "ELF3", "BIRC5"],
                  "alias_symbol": ["", "BETA-3|EEF-1", "", "", "EPR-1X", "EPR-1X"],
                  "prev_symbol": ["", "EEF1B|EEF1B3", "SNHG2", "", "", ""]}).to_csv(path, sep="\t", index=False)
    a = literature.hgnc_aliases(path)
    assert a["EEF1B2P2"] == ["EEF1B3"]      # EEF1B: stem of EEF1B2; BETA-3, EEF-1 generic
    assert a["GAS5"] == ["SNHG2"]           # SNHG20 is another number, not an extension of a stem
    assert a["BIRC5"] == [] and a["ELF3"] == []         # an alias of two genes


def test_stem_leaks():
    assert literature.stem_leaks("IGF2-AS", 1362, 2000)
    assert not literature.stem_leaks("HNF1A-AS1", 60, 1260)
    assert not literature.stem_leaks("FOO-AS1", 5, 8)         # too few papers on the stem to tell
    assert not literature.stem_leaks("LINC01090", 50, 50)     # not an antisense-style name


def test_d5_drops_leaky_names_and_dates_the_cache(tmp_path, monkeypatch):
    counts = {"IGF2-AS": 1300, "IGF2": 2000, "IGF2AS": 30}

    def hit_count(query):
        names = re.findall(r'TITLE:"([^"]+)"', query)
        return max(counts.get(n, 0) for n in names)
    monkeypatch.setattr(literature, "hit_count", hit_count)
    cache = tmp_path / "cache.json"
    out = literature.d5(pd.Series({"L": "IGF2-AS", "M": "MSTRG.1"}), {"IGF2-AS": ["IGF2AS"]}, cache, threads=1)
    assert out.loc["L", "names_dropped"] == "IGF2-AS" and out.loc["L", "aliases_used"] == "IGF2AS"
    assert out.loc["L", "n_liver"] == 30 and out.loc["M", "n_liver"] == 0
    saved = json.loads(cache.read_text())
    assert all(set(v) == {"n", "date"} for v in saved.values())
