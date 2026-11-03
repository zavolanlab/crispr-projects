"""
D5 prior knowledge (LOGIC.md, Step 3): Europe PMC counts of papers on a locus with liver co-terms.

Query per feature: the exact symbol or one of its usable HGNC previous symbols / aliases (hgnc_aliases), in title
or abstract, AND an RNA term (guards against symbols that are also words, e.g. NORAD), AND a liver term:
  n_any            symbol + RNA term
  n_liver          symbol + RNA term + liver term
  n_review         ... restricted to reviews
  n_retracted      ... restricted to retracted publications
  n_liver_cancer   n_liver papers with a cancer term
  n_liver_chronic  n_liver papers with a chronic liver disease term
Statistic D5 = log2(1 + n_liver - n_retracted); its parts for the summary scores D5_hcc = log2(1 + n_liver_cancer)
and D5_chronic = log2(1 + n_liver_chronic). Exact-phrase matching keeps papers on a neighbouring
symbol (LINC01089) from being credited to another (LINC01090). Loci named only by an Ensembl, StringTie or
NONCODE id (UNNAMED) get counts of 0 without a query (nobody publishes on an unnamed locus).

Europe PMC splits a phrase at hyphens, so "IGF2-AS" also matches "IGF2 as ...": a name of the form <gene>-AS/-IT/
-DT[n] is dropped from the query when it has at least STEM_LEAK times the liver papers of <gene> itself
(stem_leaks); a locus with no name left gets no counts (missing, not 0). The dropped names are reported.
"""
import datetime
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
# aliases that are generic tokens rather than names: a Greek letter or a short word with a number (BETA-3, EEF-1)
GENERIC_ALIAS = re.compile(r"^((alpha|beta|gamma|delta|epsilon|kappa|lambda|sigma|theta|omega)\W?\d+\w?"
                           r"|[A-Za-z]{1,3}-\d+)$", re.IGNORECASE)
LEAKY_NAME = re.compile(r"^(?P<stem>[A-Za-z0-9.]+)-(AS|IT|DT)\d*$")
STEM_LEAK = 0.5         # name dropped if it has >= this fraction of the liver papers of its stem gene ...
STEM_MIN_PAPERS = 20    # ... and the stem gene has at least this many
COUNTS = ("n_any", "n_liver", "n_review", "n_retracted", "n_liver_cancer", "n_liver_chronic")    # build_queries


def usable_alias(alias, approved):
    """Aliases that are words or short gene-like tokens (ALIEN, LEU2) collect unrelated papers: keep aliases of
    >= 5 characters that contain a digit and are not the approved symbol of another gene or a generic token
    (GENERIC_ALIAS)."""
    return (len(alias) >= 5 and any(ch.isdigit() for ch in alias) and alias not in approved
            and not GENERIC_ALIAS.match(alias))


def hgnc_aliases(path):
    """
    symbol -> list of usable aliases and previous symbols (HGNC complete set). Also left out: names that HGNC lists
    for more than one gene, and names that are the stem of another gene's symbol (EEF1B of a pseudogene, the
    common name of the EEF1B2 protein): papers using them are mostly about the other gene.
    """
    h = pd.read_csv(path, sep="\t", usecols=["symbol", "alias_symbol", "prev_symbol"], dtype=str).fillna("")
    approved = set(h["symbol"])
    names = {sym: [a for a in (prev + "|" + alias).split("|") if a and a != sym]
             for sym, alias, prev in h.itertuples(index=False)}
    owners = {}
    for sym, ns in names.items():
        for a in ns:
            owners.setdefault(a, set()).add(sym)
    by_stem = {}       # stem ending in a letter -> approved symbols that extend it by 1-2 characters (EEF1B -> EEF1B2;
    for sym in approved:                             # not SNHG2 -> SNHG20, a different number)
        for k in (1, 2):
            if sym[:-k][-1:].isalpha():
                by_stem.setdefault(sym[:-k], set()).add(sym)
    return {sym: [a for a in ns if usable_alias(a, approved) and owners[a] == {sym}
                  and not (by_stem.get(a, set()) - {sym})]
            for sym, ns in names.items()}


def stem_leaks(name, n_liver_name, n_liver_stem):
    """True if a <gene>-AS/-IT/-DT name collects the papers of <gene> (module docstring)."""
    return (LEAKY_NAME.match(name) is not None and n_liver_stem >= STEM_MIN_PAPERS
            and n_liver_name >= STEM_LEAK * n_liver_stem)


def _field(term):
    return f'(TITLE:"{term}" OR ABSTRACT:"{term}")'


def build_queries(symbol, aliases=()):
    """Europe PMC query string of each of COUNTS."""
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
    """Number of Europe PMC hits of a query (None if every attempt fails)."""
    url = EUROPEPMC + "?" + urllib.parse.urlencode({"query": query, "format": "json", "pageSize": 1,
                                                     "resultType": "idlist"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return int(json.load(r)["hitCount"])
        except Exception:
            time.sleep(2 ** attempt)
    return None


def _read_cache(path):
    """query -> {"n": hits, "date": query date}; counts cached before dates were kept have date None."""
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        cache = json.load(fh)
    return {q: v if isinstance(v, dict) else {"n": v, "date": None} for q, v in cache.items()}


def _write_cache(cache, path):
    tmp = f"{path}.tmp.{os.getpid()}"         # a kill while writing must not corrupt the cache
    with open(tmp, "w") as fh:
        json.dump(cache, fh)
    os.replace(tmp, path)


def _count(queries, cache, cache_path, threads):
    """Hit counts of the queries not yet in the cache, added to it (and saved every 200 queries)."""
    todo = sorted(set(queries) - set(cache))
    today = str(datetime.date.today())
    with ThreadPoolExecutor(threads) as ex:
        for i, (q, n) in enumerate(zip(todo, ex.map(hit_count, todo))):
            if n is not None:
                cache[q] = {"n": n, "date": today}
            if i % 200 == 199:
                _write_cache(cache, cache_path)
    _write_cache(cache, cache_path)
    failed = [q for q in todo if q not in cache]
    if failed:     # a missing count would score as 'no data' (0.5) instead of failing; rerun to retry these only
        raise RuntimeError(f"{len(failed)} Europe PMC queries failed after retries (cached counts kept), e.g. {failed[0]}")


def d5(names, aliases, cache_path, threads=4, max_aliases=5):
    """
    COUNTS and D5 per feature. `names`: Series feature_id -> symbol; `aliases`: hgnc_aliases. Counts are cached per
    query string in `cache_path` (JSON, with the query date) so that reruns only query what is new.
    """
    cache = _read_cache(cache_path)
    named = {fid: [sym, *aliases.get(sym, [])[:max_aliases]] for fid, sym in names.items()
             if isinstance(sym, str) and not UNNAMED.match(sym)}

    # names of the form <gene>-AS... against <gene> itself
    leak_checks = {n: LEAKY_NAME.match(n).group("stem") for ns in named.values() for n in ns if LEAKY_NAME.match(n)}
    liver = {n: build_queries(n)["n_liver"] for n in set(leak_checks) | set(leak_checks.values())}
    _count(liver.values(), cache, cache_path, threads)
    leaky = {n for n, stem in leak_checks.items()
             if stem_leaks(n, cache[liver[n]]["n"], cache[liver[stem]]["n"])}

    queries = {}
    for fid, ns in named.items():
        kept = [n for n in ns if n not in leaky]
        if kept:                                # else no counts: every name collects another gene's papers
            queries[fid] = build_queries(kept[0], kept[1:])
    _count({q for qs in queries.values() for q in qs.values()}, cache, cache_path, threads)

    rows = []
    for fid, sym in names.items():
        rec = {"feature_id": fid, "symbol": sym}
        if fid in named:
            rec["aliases_used"] = ",".join(n for n in named[fid][1:] if n not in leaky)
            rec["names_dropped"] = ",".join(n for n in named[fid] if n in leaky)
            for k, q in queries.get(fid, {}).items():
                rec[k] = cache[q]["n"]
        else:
            rec.update(dict.fromkeys(COUNTS, 0), aliases_used="", names_dropped="")
        rows.append(rec)
    out = pd.DataFrame(rows).set_index("feature_id").reindex(columns=["symbol", "aliases_used", "names_dropped",
                                                                      *COUNTS])
    for k in COUNTS:
        out[k] = pd.to_numeric(out[k], errors="coerce")
    out["D5"] = np.log2(1 + (out["n_liver"] - out["n_retracted"]).clip(lower=0))
    out["D5_hcc"] = np.log2(1 + out["n_liver_cancer"])
    out["D5_chronic"] = np.log2(1 + out["n_liver_chronic"])
    return out
