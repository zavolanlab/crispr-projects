"""
D3 somatic copy number (LOGIC.md, Step 3), TCGA-LIHC primary tumours.

Masked copy-number segments (GDC, GRCh38; log2 copy ratio). Per feature and tumour, the segment value at the
midpoint of the feature span. Statistics: Spearman rho of CN vs expression, partial on tumour purity (ABSOLUTE,
TCGA PanCanAtlas: stromal and immune cells shift expression and CN amplitude together; the unadjusted rho is
reported as cn_rho_raw), and gain / loss frequency (segment value > +0.3 / < -0.3).

Dosage response relative to the neighbourhood: in a frequently gained or lost region most expressed genes follow
the copy number to some degree, so cn_rho alone mostly says where a locus lies. The same purity-adjusted rho is
computed for the protein-coding genes (recount3 gene sums of the same tumours, scaled by the size factors) and
cn_rho_neighbours is the median over the NEIGHBOURS nearest of them on the same chromosome (detected, ncov > 0, in
>= MIN_DETECTED of the tumours); cn_rho_excess = cn_rho - cn_rho_neighbours.

D3 = max(cn_rho_excess, 0) x (gain frequency - loss frequency), and 0 where cn_fdr >= 0.05: > 0 when expression
follows a frequent copy gain more closely than the genes around it (oncogene-like), < 0 when it follows a frequent
loss.
"""
import numpy as np
import pandas as pd
from scipy import stats

from .disease import bh

CN_THRESHOLD = 0.3
NEIGHBOURS = 20
MIN_DETECTED = 0.2


def read_purity(path):
    """ABSOLUTE tumour purity (TCGA PanCanAtlas mastercalls) by 16-character sample barcode."""
    p = pd.read_csv(path, sep="\t")
    p = p[p["call status"] == "called"]
    return p.assign(b=p["sample"].str[:16]).drop_duplicates("b").set_index("b")["purity"].astype(float)


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


def gene_spans(exons, gene_ids):
    """chrom / start0 / end1 of the genes `gene_ids` (unversioned) from GTF exons (quant.read_gtf_exons)."""
    e = exons.assign(gene=exons["gene_id"].str.split(".").str[0])
    e = e[e["gene"].isin(set(gene_ids))]
    return e.groupby("gene").agg(chrom=("chrom", "first"), start0=("start0", "min"), end1=("end1", "max"))


def dosage_rho(expr, ncov, samples, cn, purity=None, min_detected=MIN_DETECTED):
    """Purity-adjusted Spearman rho of copy number vs expression per gene (rows of `expr`, log scale; `ncov` linear,
    for detection) in the primary tumours; genes detected in < `min_detected` of the tumours are left out."""
    tum = tcga_expression_by_barcode(expr, samples, "01")
    detected = (tcga_expression_by_barcode(ncov, samples, "01") > 0).mean(axis=1) >= min_detected
    return _spearman_rows(tum[detected], cn.reindex(index=tum.index[detected]), purity)["rho"]


def neighbour_median(spans, ref_spans, ref_values, k=NEIGHBOURS):
    """Per row of `spans`: median of `ref_values` over the k rows of `ref_spans` (same chromosome) with a value
    whose midpoints are nearest to its midpoint; NaN on chromosomes with fewer than k."""
    ref = ref_spans.assign(v=ref_values.reindex(ref_spans.index)).dropna(subset=["v"])
    ref_mid = (ref["start0"] + ref["end1"]) / 2
    mid = (spans["start0"] + spans["end1"]) / 2
    out = pd.Series(np.nan, index=spans.index)
    for chrom, idx in spans.groupby("chrom").groups.items():
        on = ref["chrom"] == chrom
        if on.sum() < k:
            continue
        rm, rv = ref_mid[on].to_numpy(), ref.loc[on, "v"].to_numpy()
        for fid in idx:
            nearest = np.argpartition(np.abs(rm - mid[fid]), k - 1)[:k]
            out[fid] = np.median(rv[nearest])
    return out


def d3(expr, samples, cn, purity=None, threshold=CN_THRESHOLD, neighbour_rho=None):
    """Per feature: CN statistics and the signed D3 statistic (rho adjusted for `purity` if given, a Series by
    16-character sample barcode). `neighbour_rho`: cn_rho_neighbours per feature (neighbour_median of dosage_rho);
    without it the excess is cn_rho itself."""
    tum = tcga_expression_by_barcode(expr, samples, "01")
    out = pd.DataFrame(index=expr.index)
    cnf = cn.reindex(index=expr.index)
    out["cn_rho_raw"] = _spearman_rows(tum, cnf)["rho"]
    c = _spearman_rows(tum, cnf, purity)
    out["cn_rho"], out["cn_p"], out["cn_n"] = c["rho"], c["p"], c["n"]
    out["cn_fdr"] = bh(out["cn_p"])
    out["cn_rho_neighbours"] = (neighbour_rho.reindex(out.index) if neighbour_rho is not None
                                else pd.Series(0.0, index=out.index))
    out["cn_rho_excess"] = out["cn_rho"] - out["cn_rho_neighbours"]
    tcn = cn.reindex(index=expr.index, columns=[s for s in cn.columns if s[13:15] == "01"])
    out["cn_gain_freq"] = (tcn > threshold).sum(axis=1) / tcn.notna().sum(axis=1)
    out["cn_loss_freq"] = (tcn < -threshold).sum(axis=1) / tcn.notna().sum(axis=1)
    d3 = out["cn_rho_excess"].clip(lower=0) * (out["cn_gain_freq"] - out["cn_loss_freq"])
    out["D3"] = d3.where(out["cn_fdr"] < 0.05, 0.0).where(out["cn_rho"].notna())
    return out
