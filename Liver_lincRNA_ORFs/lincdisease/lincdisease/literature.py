"""
D5 prior knowledge (LOGIC.md, Step 3): Europe PMC counts of papers on a locus with liver co-terms.

Query per feature: the exact symbol and its usable HGNC previous symbols / aliases (usable_alias), each in
title or abstract,
AND an RNA term (guards against symbols that are also words, e.g. NORAD), AND a liver term:
  n_any        symbol + RNA term
  n_liver      symbol + RNA term + liver term
  n_review     ... restricted to reviews
  n_retracted  ... restricted to retracted publications
Statistic D5 = log2(1 + n_liver - n_retracted). Exact-phrase matching keeps papers on a neighbouring
symbol (LINC01089) from being credited to another (LINC01090). Loci named only by an Ensembl or StringTie id
get counts of 0 without a query (nobody publishes on an unnamed locus).
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

RNA_TERMS = ['lncRNA', '"long non-coding"', '"long noncoding"', '"non-coding RNA"', '"noncoding RNA"', 'RNA',
             'transcript', 'micropeptide', 'peptide', 'gene', 'expression']
LIVER_TERMS = ['liver', 'hepatic', 'hepatocellular', 'hepatocyte*', 'HCC', 'hepatoma', 'cirrhosis',
               '"liver fibrosis"', 'NAFLD', 'MASLD', 'NASH', 'MASH', 'steatosis', 'hepatitis', 'HepG2', 'Huh7']
# split of the liver papers for the HCC and chronic liver disease scores (a paper can count for both)
CANCER_TERMS = ['cancer', 'carcinoma', 'tumour', 'tumor', 'HCC', 'hepatoma', 'hepatocellular', 'oncogen*',
                'metasta*']
CHRONIC_TERMS = ['NAFLD', 'MASLD', 'NASH', 'MASH', 'steatosis', 'steatohepatitis', '"fatty liver"', 'fibrosis',
                 'cirrhosis', 'hepatitis', '"liver injury"', '"liver failure"', 'cholestasis', 'cholestatic']
UNNAMED = re.compile(r"^(ENSG\d+|MSTRG\.\d+|NONHSAG\d+)")


def usable_alias(alias, approved):
    """Aliases that are words or short gene-like tokens (ALIEN, LEU2) collect unrelated papers: keep aliases of
    >= 5 characters that contain a digit and are not the approved symbol of another gene."""
    return len(alias) >= 5 and any(ch.isdigit() for ch in alias) and alias not in approved


def hgnc_aliases(path):
    """symbol -> list of usable aliases and previous symbols (HGNC complete set)."""
    h = pd.read_csv(path, sep="\t", usecols=["symbol", "alias_symbol", "prev_symbol"], dtype=str).fillna("")
    approved = set(h["symbol"])
    out = {}
    for sym, alias, prev in h.itertuples(index=False):
        names = [a for a in (prev + "|" + alias).split("|") if a and a != sym and usable_alias(a, approved)]
        out[sym] = names
    return out


def _field(term):
    return f'(TITLE:"{term}" OR ABSTRACT:"{term}")'


def build_queries(symbol, aliases=()):
    names = " OR ".join(_field(n) for n in [symbol, *aliases])
    rna = " OR ".join(RNA_TERMS)
    liver = " OR ".join(LIVER_TERMS)
    base = f"({names}) AND ({rna})"
    liver_q = f"{base} AND ({liver})"
    return {
        "n_any": base,
        "n_liver": liver_q,
        "n_review": f'{liver_q} AND PUB_TYPE:"review"',
        "n_retracted": f'{liver_q} AND PUB_TYPE:"Retracted Publication"',
        "n_liver_cancer": f"{liver_q} AND ({' OR '.join(CANCER_TERMS)})",
        "n_liver_chronic": f"{liver_q} AND ({' OR '.join(CHRONIC_TERMS)})",
    }


def hit_count(query, retries=4):
    url = EUROPEPMC + "?" + urllib.parse.urlencode({"query": query, "format": "json", "pageSize": 1,
                                                     "resultType": "idlist"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return int(json.load(r)["hitCount"])
        except Exception:
            time.sleep(2 ** attempt)
    return None


def d5(names, aliases, cache_path, threads=4, max_aliases=5):
    """
    `names`: Series feature_id -> symbol. Counts are cached per query string in `cache_path` (JSON) so that
    reruns only query what is new.
    """
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as fh:
            cache = json.load(fh)
    queries = {}
    for fid, sym in names.items():
        if not isinstance(sym, str) or UNNAMED.match(sym):
            continue
        queries[fid] = build_queries(sym, aliases.get(sym, [])[:max_aliases])
    todo = sorted({q for qs in queries.values() for q in qs.values()} - set(cache))
    with ThreadPoolExecutor(threads) as ex:
        for i, (q, n) in enumerate(zip(todo, ex.map(hit_count, todo))):
            if n is not None:
                cache[q] = n
            if i % 200 == 199:
                with open(cache_path, "w") as fh:
                    json.dump(cache, fh)
    with open(cache_path, "w") as fh:
        json.dump(cache, fh)
    failed = [q for q in todo if q not in cache]
    if failed:     # a missing count would score as 'no data' (0.5) instead of failing; rerun to retry these only
        raise RuntimeError(f"{len(failed)} Europe PMC queries failed after retries (cached counts kept), e.g. {failed[0]}")

    rows = []
    for fid, sym in names.items():
        rec = {"feature_id": fid, "symbol": sym}
        if fid in queries:
            rec["aliases_used"] = ",".join(aliases.get(sym, [])[:max_aliases])
            for k, q in queries[fid].items():
                rec[k] = cache.get(q)
        else:
            rec.update(n_any=0, n_liver=0, n_review=0, n_retracted=0, n_liver_cancer=0, n_liver_chronic=0,
                       aliases_used="")
        rows.append(rec)
    out = pd.DataFrame(rows).set_index("feature_id")
    for k in ("n_any", "n_liver", "n_review", "n_retracted", "n_liver_cancer", "n_liver_chronic"):
        out[k] = pd.to_numeric(out[k], errors="coerce")
    out["D5"] = np.log2(1 + (out["n_liver"] - out["n_retracted"]).clip(lower=0))
    out["D5_hcc"] = np.log2(1 + out["n_liver_cancer"])
    out["D5_chronic"] = np.log2(1 + out["n_liver_chronic"])
    return out
