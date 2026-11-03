"""
D1 chronic liver disease (LOGIC.md, Step 3): expression of loci and controls in recount3 SRA cohorts of
MASLD/MASH and viral hepatitis liver biopsies, quantified like TCGA/GTEx (same regions, `ncov`).

Cohorts (phenotypes from the recount3 SRA `sample_attributes`):
  hoang   SRP197353 (GSE130970, Hoang et al. 2019): 78 NAFLD biopsies, fibrosis stage 0-4 and NAFLD activity score
  suppli  SRP186450 (GSE126848, Suppli et al. 2019): 14 normal-weight, 12 obese, 15 NAFL, 16 NASH
  hcv     SRP078450 (GSE84346, Boldanova et al. 2017): untreated biopsies, 22 HCV-infected, 6 non-HCV
  govaere GSE135251 (Govaere et al. 2020): 206 NAFLD biopsies (NAFL, NASH F0-F4) and 10 controls, fibrosis stage and
          NAS; aligned here (STAR, GENCODE v42; workflow/rules/govaere.smk), loci counted in the CRAMs over the same
          regions (samtools bedcov), gene counts and technical covariates from STAR (exonic fraction = reads
          assigned to genes / uniquely mapped reads)

Sample scaling: `ncov` is divided by median-of-ratios size factors of the recount3 gene sums, per cohort
(expression.py; hepatocyte transcripts, which dominate liver RNA, lose their share with injury).

Per cohort, a signed partial rank correlation in [-1, 1] of log2(scaled ncov + pseudocount) with the phenotype
(> 0 = higher with disease): fibrosis stage (hoang, govaere; NAS reported), NAFL + NASH vs normal-weight + obese (suppli),
HCV vs non-HCV (hcv). Ranks of expression and phenotype are residualized on the ranks of the technical covariates
(log input reads, % uniquely mapped reads, exonic fraction; Suppli cases differ from controls in all three), and
the p-value is a t test with n - 2 - k degrees of freedom.
  scored         technical covariates; for loci whose expression follows a disease-associated parent gene
                 (`parent_driven`), the parent gene's expression is added as a covariate
  cell-intrinsic technical covariates + cell-type marker scores (reported, not scored: a lincRNA of stellate
                 cells that rises with fibrosis is disease-associated, but not regulated within cells)
D1_chronic = weighted mean of the available scored cohort statistics, per feature: fibrosis 2/3, split between hoang
and govaere in proportion to the samples tested for the feature (all of it to one cohort if the other has no
statistic); hcv 1/3; suppli 0; a feature without a fibrosis statistic is the hcv statistic, and vice versa. The combined
p-value is Stouffer's (z from the two-sided p with the sign of the statistic, weights sqrt(n)), BH-FDR over features.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

from . import quant
from .disease import bh
from .expression import size_factors

COHORTS = {"hoang": "SRP197353", "suppli": "SRP186450", "hcv": "SRP078450"}    # recount3 projects
ALL_COHORTS = [*COHORTS, "govaere"]
# Cohort weights of D1_chronic (cohort_weights). The fibrosis cohorts share FIBROSIS_WEIGHT, per feature in
# proportion to the samples tested; the other cohorts have fixed weights (WEIGHTS). Suppli has weight 0 (reported,
# not scored): its cases differ technically from its controls (fewer uniquely mapped reads, lower exonic fraction;
# the 15 samples with the lowest exonic fraction are all cases), and after scaling and covariate adjustment 40% of
# loci still have p < 0.01 with a median r of +0.14 (Hoang: 9%, 0.00).
FIBROSIS_COHORTS = ("hoang", "govaere")
FIBROSIS_WEIGHT = 2 / 3
WEIGHTS = {"suppli": 0.0, "hcv": 1 / 3}
MIN_DETECTED = 0.2      # features detected (ncov > 0) in fewer samples of a cohort are not tested
TECHNICAL = ["log_input_reads", "pct_unique", "exonic_fraction"]
# liver cell-type markers (MacParland et al. 2018 liver atlas); COL1A1 left out, as it is the D1 chronic control,
# and KRT8/KRT18, as processed pseudogenes of them are among the loci
CELL_MARKERS = {
    "hepatocyte": ["CYP2E1", "CYP3A4", "CYP1A2", "APOC3", "TF", "SERPINA1", "HPX", "F9"],
    "stellate": ["DCN", "LUM", "COL3A1", "PDGFRB", "RGS5", "LRAT", "RELN"],
    "cholangiocyte": ["KRT7", "KRT19", "EPCAM", "SOX9", "CFTR"],
    "macrophage": ["CD68", "CD163", "MARCO", "VSIG4", "C1QA"],
    "lymphocyte": ["PTPRC", "CD3E", "CD2", "NKG7", "MS4A1"],
    "endothelial": ["CLEC4G", "STAB2", "FCN3", "PECAM1", "CLEC1B"],
}
# loci whose expression correlates this strongly (median Spearman over cohorts) with a parent gene that is itself
# associated with chronic disease (FDR < 0.05, same sign) are parent-adjusted; at 0.8 YWHAZP2/3 (rho 0.6-0.63 with
# YWHAZ, D1_chronic 0.66) escaped while their own signal (0.53-0.59) was mostly the parent's
PARENT_RHO = 0.5


def parse_attributes(values):
    """recount3 `sample_attributes` ('key;;value|key;;value') -> DataFrame of attributes."""
    rows = []
    for v in values:
        d = {}
        for item in (v.split("|") if isinstance(v, str) else []):
            k, _, val = item.partition(";;")
            d[k.strip()] = val.strip()
        rows.append(d)
    return pd.DataFrame(rows)


def phenotypes(man):
    """Manifest rows (one recount3 SRA project each) -> cohort, fibrosis_stage, nas, group (case/control/NaN)."""
    a = parse_attributes(man["sample_attributes"]).set_index(man.index)
    out = pd.DataFrame(index=man.index)
    out["cohort"] = man["project"].map({v: k for k, v in COHORTS.items()})
    get = lambda col: a[col] if col in a else pd.Series(np.nan, index=a.index)
    out["fibrosis_stage"] = pd.to_numeric(get("fibrosis stage"), errors="coerce")
    out["nas"] = pd.to_numeric(get("nafld activity score"), errors="coerce")
    disease = get("disease").str.lower()
    infection = get("infection status")
    treated = get("treatment").fillna("None") != "None"
    group = pd.Series(np.nan, index=man.index, dtype=object)
    s = out["cohort"] == "suppli"
    group[s & disease.isin(["nafld", "nash"])] = "case"
    group[s & disease.isin(["healthy", "obese"])] = "control"
    h = (out["cohort"] == "hcv") & ~treated
    group[h & (infection == "HCV")] = "case"
    group[h & (infection == "non-HCV")] = "control"
    out["group"] = group
    return out


def manifest():
    """Samples of the chronic-disease cohorts (quant.manifest format) with their phenotypes."""
    parts = [quant.manifest("sra", p, ["sample_attributes", "sample_title"]) for p in COHORTS.values()]
    man = pd.concat(parts, ignore_index=True)
    man = pd.concat([man, phenotypes(man)], axis=1)
    keep = man["fibrosis_stage"].notna() | man["group"].notna()     # drops IFN-treated HCV biopsies
    return man[keep].reset_index(drop=True)



def read_soft(soft):
    """GEO sample annotation (acc.cgi text format) -> one row per GSM with its sample characteristics as columns."""
    rec, cur = [], None
    for line in open(soft):
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=")[1].strip()}
            rec.append(cur)
        elif line.startswith("!Sample_characteristics_ch1") and cur is not None:
            k, _, v = line.split("=", 1)[1].strip().partition(": ")
            cur[k] = v
    return pd.DataFrame(rec)


def govaere_samples(runs_tsv, soft):
    """Govaere runs (index run accession) with their GEO sample characteristics."""
    runs = pd.read_csv(runs_tsv, sep="\t")[["run_accession", "experiment_alias"]]
    return runs.merge(read_soft(soft), left_on="experiment_alias", right_on="gsm").set_index("run_accession")


def govaere(gov_dir, runs_tsv, soft, regions, exclude=()):
    """
    Govaere runs of the run table except `exclude` -> ncov (features x samples), STAR gene counts (unversioned
    gene id x samples; column 4, reverse-stranded library), QC rows in recount_qc column names, phenotypes.
    Sample ids 'govaere:<run>'. `soft`: GEO sample annotation with the sample characteristics. Raises if a run has
    no locus coverage or the coverage was computed on regions other than `regions` (current quant/regions.tsv).
    """
    from .expression import changed_features
    meta = govaere_samples(runs_tsv, soft)
    done = [r for r in meta.index if r not in set(exclude)]
    cov_dir = os.path.join(gov_dir, "locus_cov")
    missing = [r for r in done if not os.path.exists(os.path.join(cov_dir, f"{r}.tsv"))]
    if missing:
        raise ValueError(f"no locus coverage for {len(missing)} Govaere runs ({missing[:3]})")
    saved = os.path.join(cov_dir, "regions.tsv")
    if not os.path.exists(saved) or changed_features(pd.read_csv(saved, sep="\t"), regions):
        raise ValueError(f"{cov_dir} was computed on other regions than quant/regions.tsv; delete it to recompute")
    ids = ["govaere:" + r for r in done]
    # summed per-base depth over the regions of each feature: proportional to ncov within the cohort, and the size
    # factors are relative within each cohort (d1_chronic)
    ncov = pd.concat([pd.read_csv(os.path.join(cov_dir, f"{r}.tsv"), sep="\t", index_col=0).iloc[:, 0]
                      for r in done], axis=1, keys=ids)
    # a feature is absent from a run's file only if the file predates it (bedcov reports every region)
    stale = ncov.columns[ncov.isna().any()]
    if len(stale):
        raise ValueError(f"{len(stale)} Govaere locus_cov files lack features of the others (older regions.tsv); "
                         f"delete and recompute them: {', '.join(c.split(':')[1] for c in stale[:10])}")
    genes = pd.concat([pd.read_csv(os.path.join(gov_dir, "counts", f"{r}.ReadsPerGene.tsv"), sep="\t", header=None,
                                   index_col=0)[3] for r in done], axis=1, keys=ids)
    assigned = genes.loc[~genes.index.str.startswith("N_")].sum()
    genes = genes.loc[~genes.index.str.startswith("N_")]
    genes.index = genes.index.str.split(".").str[0]
    genes = genes.groupby(level=0).sum()
    qc = {}
    for r, i in zip(done, ids):
        d = dict((k.strip(), v.strip()) for k, v in (line.split("|", 1) for line in
                 open(os.path.join(gov_dir, "star_logs", f"{r}.Log.final.out")) if "|" in line))
        qc[i] = {"star.number_of_input_reads_both": int(d["Number of input reads"]),
                 "star.uniquely_mapped_reads_%_both": float(d["Uniquely mapped reads %"].rstrip("%")),
                 "exonic_fraction": assigned[i] / int(d["Uniquely mapped reads number"])}
    pheno = pd.DataFrame({"cohort": "govaere",
                          "fibrosis_stage": pd.to_numeric(meta.loc[done, "fibrosis stage"], errors="coerce").to_numpy(),
                          "nas": pd.to_numeric(meta.loc[done, "nas score"], errors="coerce").to_numpy(),
                          "group": np.nan}, index=ids)
    return ncov, genes, pd.DataFrame(qc).T, pheno


def technical_covariates(qc, genes_all):
    """
    recount_qc rows (index sample_id) and all gene sums -> TECHNICAL covariates per sample. Exonic fraction: from
    `qc` where given (Govaere: STAR reads assigned to genes / uniquely mapped reads), otherwise the sum of the
    AUC-scaled recount3 gene sums / 1e9 (fraction of the base coverage in genes).
    """
    ef = genes_all.sum() / 1e9
    if "exonic_fraction" in qc:
        ef = pd.to_numeric(qc["exonic_fraction"], errors="coerce").reindex(ef.index).fillna(ef)
    return pd.DataFrame({
        "log_input_reads": np.log10(pd.to_numeric(qc["star.number_of_input_reads_both"], errors="coerce")),
        "pct_unique": pd.to_numeric(qc["star.uniquely_mapped_reads_%_both"], errors="coerce"),
        "exonic_fraction": ef,
    })


def cell_scores(log_genes, symbol_to_id):
    """Mean z-score of the log expression of each marker set (sets with < 3 expressed markers are left out)."""
    out = {}
    for cell, genes in CELL_MARKERS.items():
        ids = [symbol_to_id[g] for g in genes if g in symbol_to_id and symbol_to_id[g] in log_genes.index]
        if len(ids) < 3:
            continue
        x = log_genes.loc[ids]
        out[cell] = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1).replace(0, np.nan), axis=0).mean()
    return pd.DataFrame(out)


def partial_rank(expr, y, covars=None):
    """
    Partial rank correlation of every row of `expr` (features x samples) with y, given the sample covariates
    (samples x k; ranks are residualized on covariate ranks). Samples with a missing y or covariate are dropped.
    Returns DataFrame stat / p / n.
    """
    ok = y.notna()
    if covars is not None and covars.shape[1]:
        ok &= covars.notna().all(axis=1).reindex(y.index, fill_value=False)
    idx = y.index[ok]
    n = len(idx)
    k = 0 if covars is None else covars.shape[1]
    if n - 2 - k < 3:
        return pd.DataFrame(columns=["stat", "p", "n"], dtype=float)
    rx = expr[idx].rank(axis=1).to_numpy(float)
    ry = y[idx].rank().to_numpy(float)
    X = np.ones((n, 1))
    if k:
        X = np.column_stack([X, covars.loc[idx].rank().to_numpy(float)])
    H = X @ np.linalg.pinv(X)
    ex = rx - rx @ H.T
    ey = ry - H @ ry
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (ex @ ey) / np.sqrt((ex ** 2).sum(axis=1) * (ey ** 2).sum())
        df = n - 2 - k
        t = r * np.sqrt(df / np.clip(1 - r ** 2, 1e-12, None))
    p = 2 * stats.t.sf(np.abs(t), df)
    return pd.DataFrame({"stat": r, "p": p, "n": float(n)}, index=expr.index)


def _phenotype(pheno, idx, cohort):
    if cohort in FIBROSIS_COHORTS:
        return pheno.loc[idx, "fibrosis_stage"].astype(float)
    g = pheno.loc[idx, "group"]
    return (g == "case").astype(float).where(g.notna())


def cohort_stats(expr, ncov, pheno, covars, suffix=""):
    """Per feature and cohort: partial rank statistic, p, n (columns <cohort>_r/_p/_n<suffix>); hoang also vs NAS.
    covars: {cohort: samples x k DataFrame}."""
    res = {}
    for cohort in ALL_COHORTS:
        idx = pheno.index[pheno["cohort"] == cohort]
        if not len(idx):
            continue
        detected = (ncov[idx] > 0).mean(axis=1) >= MIN_DETECTED
        e = expr.loc[detected, idx]
        z = covars.get(cohort)
        t = partial_rank(e, _phenotype(pheno, idx, cohort), z)
        if cohort in FIBROSIS_COHORTS:
            nas = partial_rank(e, pheno.loc[idx, "nas"].astype(float), z)
            res[f"{cohort}_nas_rho{suffix}"], res[f"{cohort}_nas_p{suffix}"] = nas["stat"], nas["p"]
        res[f"{cohort}_r{suffix}"], res[f"{cohort}_p{suffix}"], res[f"{cohort}_n{suffix}"] = t["stat"], t["p"], t["n"]
    return pd.DataFrame(res).reindex(expr.index)


def cohort_weights(cs, weights=WEIGHTS, suffix=""):
    """Per feature and scored cohort: weight of its statistic (missing statistic -> 0). The fibrosis cohorts share
    FIBROSIS_WEIGHT in proportion to the samples tested for the feature."""
    w = pd.DataFrame(index=cs.index)
    fib = [c for c in FIBROSIS_COHORTS if f"{c}_r{suffix}" in cs]
    n = pd.DataFrame({c: cs[f"{c}_n{suffix}"].where(np.isfinite(cs[f"{c}_r{suffix}"]), 0).fillna(0) for c in fib},
                     index=cs.index)
    tot = n.sum(axis=1)
    for c in fib:
        w[c] = FIBROSIS_WEIGHT * (n[c] / tot.where(tot > 0)).fillna(0)
    for c, wc in weights.items():
        if wc > 0 and f"{c}_r{suffix}" in cs:
            w[c] = np.where(np.isfinite(cs[f"{c}_r{suffix}"]), wc, 0.0)
    return w


def combine(cs, weights=WEIGHTS, suffix=""):
    """D1_chronic (weighted mean of the cohort statistics, cohort_weights), Stouffer p and BH-FDR."""
    w = cohort_weights(cs, weights, suffix)
    cohorts = list(w.columns)
    r = cs[[f"{c}_r{suffix}" for c in cohorts]].to_numpy(float)
    p = cs[[f"{c}_p{suffix}" for c in cohorts]].to_numpy(float)
    n = cs[[f"{c}_n{suffix}" for c in cohorts]].to_numpy(float)
    w = w.to_numpy(float)
    ok = np.isfinite(r) & (w > 0)                     # cohorts that contribute to a feature
    w_ok = np.where(ok, w, 0.0)
    wsum = w_ok.sum(axis=1)
    weighted_r = (np.where(ok, r, 0.0) * w_ok).sum(axis=1)
    out = pd.DataFrame(index=cs.index)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["D1_chronic"] = np.where(wsum > 0, weighted_r / wsum, np.nan)
    z = np.sign(r) * stats.norm.isf(np.clip(p, 1e-300, 1) / 2)
    sw = np.where(ok, np.sqrt(n), 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        zc = np.nansum(np.where(ok, z * sw, 0), axis=1) / np.sqrt((sw ** 2).sum(axis=1))
    zc = np.where(ok.any(axis=1), zc, np.nan)
    out["chronic_z"] = zc
    out["chronic_p"] = np.where(np.isfinite(zc), 2 * stats.norm.sf(np.abs(zc)), np.nan)
    out["chronic_fdr"] = bh(out["chronic_p"])
    out["chronic_n_cohorts"] = ok.sum(axis=1)
    return out


def _parent_adjusted(expr, ncov, pheno, covars, parent_expr, parents, suffix=""):
    """Cohort statistics of loci with a parent gene, with the parent's expression as an extra covariate."""
    rows = []
    for fid, gid in parents.items():
        if fid not in expr.index or gid not in parent_expr.index:
            continue
        cv = {c: z.assign(parent=parent_expr.loc[gid, z.index]) for c, z in covars.items()}
        rows.append(cohort_stats(expr.loc[[fid]], ncov.loc[[fid]], pheno, cv, suffix=suffix))
    return pd.concat(rows) if rows else pd.DataFrame()


def _parent_check(expr, genes, log_genes, pheno, tech, parents, scored):
    """
    Parent genes of the loci (`parents`: locus -> gene id): parent expression, the parents' own D1_chronic, and per
    locus parent_gene / parent_rho (median over cohorts of the Spearman correlation of locus and parent) /
    parent_D1_chronic / parent_driven (rho >= PARENT_RHO, parent FDR < 0.05, same sign as the locus' `scored`
    D1_chronic). Returns (table indexed like `scored`, parent log expression, the parents of the driven loci).
    """
    parents = parents[parents.isin(genes.index) & parents.index.isin(expr.index)]
    parent_expr = log_genes.loc[parents.unique()]
    # the parents' own association; their gene sums decide detection, as ncov does for the loci
    parent_assoc = combine(cohort_stats(parent_expr, genes.loc[parents.unique()], pheno, tech))
    cohort_samples = [pheno.index[pheno["cohort"] == c] for c in ALL_COHORTS]
    rho = pd.Series({fid: np.nanmedian([stats.spearmanr(expr.loc[fid, idx], parent_expr.loc[gid, idx])[0]
                                        for idx in cohort_samples])
                     for fid, gid in parents.items()})
    out = pd.DataFrame(index=scored.index)
    out["parent_gene"] = parents.reindex(out.index)
    out["parent_rho"] = rho.reindex(out.index)
    out["parent_D1_chronic"] = parents.map(parent_assoc["D1_chronic"]).reindex(out.index)
    parent_fdr = parents.map(parent_assoc["chronic_fdr"]).reindex(out.index)
    driven = ((out["parent_rho"] >= PARENT_RHO) & (parent_fdr < 0.05)
              & (np.sign(out["parent_D1_chronic"]) == np.sign(scored["D1_chronic"])))
    out["parent_driven"] = driven.fillna(False)
    return out, parent_expr, parents[out.index[out["parent_driven"]]]


def d1_chronic(ncov, pheno, genes, qc, symbol_to_id, protein_coding, parents=None, pseudocount=0.05):
    """
    ncov: features x samples (AUC-scaled); pheno: phenotypes() indexed by sample_id; genes: AUC-scaled recount3
    gene sums (unversioned gene id x samples); qc: recount_qc rows indexed by sample_id; symbol_to_id: gene symbol ->
    unversioned id; protein_coding: ids of protein-coding genes (size factors); parents: locus -> parent gene id
    (homology), for the parent check. Returns the feature table and the size factors.

    1. per cohort: size factors, technical covariates, cell-type marker scores;
    2. scored statistics (technical covariates) and cell-intrinsic ones (+ cell-type markers), combined over cohorts;
    3. parent check: loci driven by a disease-associated parent gene get both statistics recomputed with the parent's
       expression as an extra covariate, and the combination and FDR are redone over all features.
    """
    ncov = ncov[pheno.index.intersection(ncov.columns).intersection(genes.columns)]
    pheno = pheno.loc[ncov.columns]
    genes = genes[ncov.columns]

    # 1.
    pc = genes.loc[genes.index.intersection(pd.Index(protein_coding))]
    sf, tech, cells = {}, {}, {}
    log_genes = np.log2(genes + 1)          # samples of every cohort are overwritten with scaled values below
    for cohort in ALL_COHORTS:
        idx = pheno.index[pheno["cohort"] == cohort]
        if not len(idx):
            continue
        sf[cohort] = size_factors(pc[idx])
        lg = np.log2(genes[idx].div(sf[cohort], axis=1) + 1)
        log_genes[idx] = lg
        tech[cohort] = technical_covariates(qc.loc[idx], genes[idx])
        cells[cohort] = pd.concat([tech[cohort], cell_scores(lg, symbol_to_id)], axis=1)
    sf = pd.concat(sf.values())
    expr = np.log2(ncov.div(sf[ncov.columns], axis=1) + pseudocount)

    # 2.
    cs = cohort_stats(expr, ncov, pheno, tech)
    cell = cohort_stats(expr, ncov, pheno, cells, suffix="_cell")
    out = pd.concat([combine(cs), cs], axis=1)
    cell_r = [c for c in cell.columns if c.endswith("_r_cell")]

    def put_cell(cell):
        cc = combine(cell, suffix="_cell")
        out["D1_chronic_cell"], out["chronic_cell_fdr"] = cc["D1_chronic"], cc["chronic_fdr"]
        out[cell_r] = cell[cell_r]

    put_cell(cell)
    out["parent_driven"] = False
    if parents is None or not len(parents):
        return out, sf

    # 3.
    check, parent_expr, driven = _parent_check(expr, genes, log_genes, pheno, tech, parents, out)
    out[check.columns] = check
    if len(driven):
        adj = _parent_adjusted(expr, ncov, pheno, tech, parent_expr, driven)
        out.loc[adj.index, adj.columns] = adj
        redo = combine(out[cs.columns])
        out[redo.columns] = redo      # FDR over all features, with the parent-adjusted statistics in place
        adj = _parent_adjusted(expr, ncov, pheno, cells, parent_expr, driven, suffix="_cell")
        cell.loc[adj.index, adj.columns] = adj
        put_cell(cell)
    return out, sf
