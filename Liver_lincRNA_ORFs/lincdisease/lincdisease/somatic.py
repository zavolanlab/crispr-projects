"""
D3 somatic genome (LOGIC.md, Step 3), TCGA-LIHC primary tumours.

Copy number: masked copy-number segments (GDC, GRCh38; log2 copy ratio). Per feature and tumour, the segment
value at the midpoint of the feature span. Statistics: Spearman rho of CN vs expression, gain / loss frequency
(segment value > +0.3 / < -0.3).
Methylation: Illumina 450K beta values (GDC, GRCh38 probe coordinates). Promoter probes lie within -1500..+500 nt
of the TSS of the dominant isoform (lincRNA loci: highest TPM in HepG2/Huh7; GENCODE controls: Ensembl canonical
transcript); the promoter beta of a sample is their mean. Statistics: Spearman rho of beta vs expression in
tumours, and median beta in adjacent normals minus median beta in tumours (delta_beta).
Both rho are partial Spearman correlations adjusted for tumour purity (ABSOLUTE, TCGA PanCanAtlas): stromal and
immune cells shift expression, CN amplitude and methylation together. The unadjusted rho are reported as *_rho_raw.

Signed parts (> 0 = expression raised by a somatic change, i.e. oncogene-like):
  d3_cn   = max(rho_CN, 0) x (gain frequency - loss frequency)
  d3_meth = max(-rho_meth, 0) x clip(delta_beta / 0.2, -1, 1)
D3 = the part with the larger magnitude.
"""
import gzip
import re

import numpy as np
import pandas as pd
from scipy import stats

from .disease import bh

CN_THRESHOLD = 0.3
PROMOTER = (1500, 500)
DELTA_BETA_SCALE = 0.2


def tss_windows(transcripts, upstream=PROMOTER[0], downstream=PROMOTER[1]):
    """`transcripts`: feature_id/chrom/start0/end1/strand per isoform -> one promoter window per isoform."""
    t = transcripts
    tss = np.where(t["strand"] == "+", t["start0"], t["end1"] - 1)
    lo = np.where(t["strand"] == "+", tss - upstream, tss - downstream)
    hi = np.where(t["strand"] == "+", tss + downstream, tss + upstream)
    return pd.DataFrame({"feature_id": t["feature_id"].to_numpy(), "chrom": t["chrom"].to_numpy(),
                         "start0": lo, "end1": hi + 1})


def canonical_transcripts(gtf_path):
    """Ids of the transcripts tagged Ensembl_canonical in a GENCODE GTF."""
    ids = set()
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if "\ttranscript\t" in line and "Ensembl_canonical" in line:
                ids.add(re.search('transcript_id "([^"]+)"', line).group(1))
    return ids


def dominant_isoforms(transcripts):
    """One isoform per feature: the canonical one if flagged (column `canonical`), else the highest `TPM_max`."""
    t = transcripts.assign(_c=transcripts["canonical"].fillna(False).astype(bool),
                           _e=transcripts["TPM_max"].astype(float).fillna(-1))
    t = t.sort_values(["feature_id", "_c", "_e"], ascending=[True, False, False])
    return t.drop_duplicates("feature_id").drop(columns=["_c", "_e"])


def promoter_probes(windows, probemap):
    """Dict feature_id -> sorted list of probe ids within any of its promoter windows."""
    pm = probemap.dropna(subset=["chrom", "chromStart"])
    by = {c: g.sort_values("chromStart") for c, g in pm.groupby("chrom")}
    out = {}
    for fid, g in windows.groupby("feature_id"):
        ids = set()
        for _, w in g.iterrows():
            p = by.get(w["chrom"])
            if p is None:
                continue
            pos = p["chromStart"].to_numpy()
            ids.update(p["id"].iloc[np.searchsorted(pos, w["start0"]):np.searchsorted(pos, w["end1"])])
        if ids:
            out[fid] = sorted(ids)
    return out


def read_probemap(path):
    pm = pd.read_csv(path, sep="\t")
    return pm.rename(columns={"#id": "id"})


def read_methylation(path, probes):
    """Rows of the probes x samples beta matrix (Xena layout) for the requested probes only."""
    want = set(probes)
    rows, ids = [], []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            pid, _, rest = line.partition("\t")
            if pid in want:
                ids.append(pid)
                rows.append([float(x) if x not in ("", "NA") else np.nan for x in rest.rstrip("\n").split("\t")])
    return pd.DataFrame(rows, index=ids, columns=header)


def read_purity(path):
    """ABSOLUTE tumour purity (TCGA PanCanAtlas mastercalls) by 16-character sample barcode."""
    p = pd.read_csv(path, sep="\t")
    p = p[p["call status"] == "called"]
    return p.assign(b=p["sample"].str[:16]).drop_duplicates("b").set_index("b")["purity"].astype(float)


def promoter_beta(beta, probes_by_feature):
    """features x samples mean promoter beta."""
    return pd.DataFrame({fid: beta.reindex(p).mean(axis=0) for fid, p in probes_by_feature.items()}).T


def cn_at(spans, segments):
    """features x samples segment value at the midpoint of each feature span (NaN outside segments)."""
    seg = segments.rename(columns={"Chrom": "chrom", "Start": "start", "End": "end"})
    seg["chrom"] = "chr" + seg["chrom"].astype(str).str.replace("^chr", "", regex=True)
    mid = ((spans["start0"] + spans["end1"]) // 2).to_numpy()
    out = {}
    for sample, s in seg.groupby("sample"):
        vals = np.full(len(spans), np.nan)
        for chrom, g in s.groupby("chrom"):
            idx = np.flatnonzero(spans["chrom"].to_numpy() == chrom)
            if not len(idx):
                continue
            g = g.sort_values("start")
            st, en, v = g["start"].to_numpy(), g["end"].to_numpy(), g["value"].to_numpy()
            k = np.searchsorted(st, mid[idx], "right") - 1
            ok = (k >= 0) & (mid[idx] <= en[np.clip(k, 0, None)])
            vals[idx[ok]] = v[k[ok]]
        out[sample] = vals
    return pd.DataFrame(out, index=spans.index)


def tcga_expression_by_barcode(expr, samples, code):
    """Expression columns of TCGA samples with sample-type `code` ('01' tumour, '11' normal), relabelled by the
    16-character sample barcode; aliquots of one sample are averaged."""
    t = samples[(samples["source"] == "tcga") & (samples["tcga_barcode"].str[13:15] == code)]
    lab = t.set_index("sample_id")["tcga_barcode"].str[:16]
    e = expr[lab.index]
    return e.T.groupby(lab.to_numpy()).mean().T


def partial_spearman(x, y, z):
    """Spearman rho of x and y adjusted for z (Pearson correlation of rank residuals), with its p (n - 3 df)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    res = lambda r: r - np.polyval(np.polyfit(rz, r, 1), rz)
    r = float(np.corrcoef(res(rx), res(ry))[0, 1])
    df = len(x) - 3
    t = r * np.sqrt(df / max(1 - r * r, 1e-12))
    return r, float(2 * stats.t.sf(abs(t), df))


def _spearman_rows(a, b, covariate=None, min_n=20):
    """Row-wise Spearman rho and p between two features x samples frames (same index); partial on `covariate`
    (a Series over the columns) if given, using the samples where it is known."""
    cols = a.columns.intersection(b.columns)
    if covariate is not None:
        cols = cols.intersection(covariate.dropna().index)
        z = covariate[cols].to_numpy(float)
    rho, p, n = {}, {}, {}
    for fid in a.index.intersection(b.index):
        x, y = a.loc[fid, cols].to_numpy(float), b.loc[fid, cols].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        n[fid] = ok.sum()
        if ok.sum() < min_n or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
            continue
        if covariate is None:
            rho[fid], p[fid] = stats.spearmanr(x[ok], y[ok])
        else:
            rho[fid], p[fid] = partial_spearman(x[ok], y[ok], z[ok])
    return pd.DataFrame({"rho": pd.Series(rho, dtype=float), "p": pd.Series(p, dtype=float),
                         "n": pd.Series(n, dtype=float)})


def d3(expr, samples, cn, beta, purity=None, threshold=CN_THRESHOLD, delta_scale=DELTA_BETA_SCALE):
    """Per feature: CN and methylation statistics and the signed D3 statistic (rho adjusted for `purity` if
    given, a Series by 16-character sample barcode)."""
    tum = tcga_expression_by_barcode(expr, samples, "01")
    out = pd.DataFrame(index=expr.index)
    cnf = cn.reindex(index=expr.index)
    out["cn_rho_raw"] = _spearman_rows(tum, cnf)["rho"]
    c = _spearman_rows(tum, cnf, purity)
    out["cn_rho"], out["cn_p"], out["cn_n"] = c["rho"], c["p"], c["n"]
    out["cn_fdr"] = bh(out["cn_p"])
    tcn = cn.reindex(index=expr.index, columns=[s for s in cn.columns if s[13:15] == "01"])
    out["cn_gain_freq"] = (tcn > threshold).sum(axis=1) / tcn.notna().sum(axis=1)
    out["cn_loss_freq"] = (tcn < -threshold).sum(axis=1) / tcn.notna().sum(axis=1)

    b = beta.reindex(index=expr.index)
    out["meth_rho_raw"] = _spearman_rows(tum, b)["rho"]
    m = _spearman_rows(tum, b, purity)
    out["meth_rho"], out["meth_p"], out["meth_n"] = m["rho"], m["p"], m["n"]
    out["meth_fdr"] = bh(out["meth_p"])
    bt = b[[s for s in b.columns if s[13:15] == "01"]].median(axis=1)
    bn = b[[s for s in b.columns if s[13:15] == "11"]].median(axis=1)
    out["meth_tumour_median"], out["meth_normal_median"] = bt, bn
    out["delta_beta"] = bn - bt

    out["d3_cn"] = out["cn_rho"].clip(lower=0) * (out["cn_gain_freq"] - out["cn_loss_freq"])
    out["d3_meth"] = (-out["meth_rho"]).clip(lower=0) * (out["delta_beta"] / delta_scale).clip(-1, 1)
    parts = out[["d3_cn", "d3_meth"]]
    pick = parts.abs().fillna(-1).to_numpy().argmax(axis=1)
    out["D3"] = np.where(parts.notna().any(axis=1), parts.to_numpy()[np.arange(len(parts)), pick], np.nan)
    out["D3_source"] = np.where(out["D3"].notna(), np.array(["cn", "meth"])[pick], "")
    return out
