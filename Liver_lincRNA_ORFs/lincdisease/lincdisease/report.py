"""
Final tables (LOGIC.md, "Outputs") and the control panel check (Step 5).

Per ORF: QC flags and tier, copy-attribution confidence, the locus' disease components and score, the
actionability score, mechanism class scores, shares and dominant class.
Actionability (per ORF) = equal-weight composite of the percentiles (among all ORFs) of the locus' tumour
specificity, the number of lincORF-specific HLA-I peptides and predicted HLA-I binding (missing until MHCflurry
predictions are added).
"""
import numpy as np
import pandas as pd

from . import scoring

ACTIONABILITY = ["tumour_specificity", "n_hla_specific", "hla_binding"]


def actionability(orfs, disease_loci, peptides):
    f = pd.DataFrame(index=orfs.index)
    f["tumour_specificity"] = orfs["locus_id"].map(disease_loci["tumour_specificity"]).to_numpy()
    f["n_hla_specific"] = orfs["orf_id"].map(peptides["n_hla_specific"]).fillna(0).to_numpy()
    f["hla_binding"] = np.nan
    pct = scoring.component_percentiles(f, f, ACTIONABILITY)
    res = scoring.composite(pct, dict.fromkeys(ACTIONABILITY, 1.0))
    return pd.concat([f.add_prefix("act_"), pct.add_prefix("act_pct_"),
                      res.rename(columns={"score": "actionability_score", "coverage": "actionability_coverage"})],
                     axis=1)


def orf_table(orfs, qc_table, disease_loci, mech, act):
    keep_disease = ["disease_score", "disease_coverage", "disease_score_noD5", "disease_direction",
                    "comp_D1", "comp_D2", "comp_D3", "comp_D4", "comp_D5", "comp_D6",
                    "log2FC_paired", "fdr_paired", "auc_vs_gtex", "fdr_vs_gtex", "HR_OS", "fdr_OS", "HR_PFI",
                    "fdr_PFI", "gwas_top_trait", "gwas_top_snp", "gwas_disease_max", "gwas_enzyme_max",
                    "n_liver", "n_review", "n_retracted", "tumour_specificity"]
    keep_disease += [c for c in ("disease_score_fitted",) if c in disease_loci]
    dl = disease_loci[keep_disease].reindex(orfs["locus_id"])
    dl.index = orfs.index
    mcols = [c for c in mech.columns if c.endswith(("_score", "_share", "_coverage")) or c == "dominant_class"]
    qcols = ["confidence_tier", "copy_attribution", "segdup_fracMatch", "grch38_false_dup", "antisense_to_coding",
             "small_RNA_overlap", "multimapper_dominated", "coding_gene_extension", "tcga_quantifiable"]
    base = orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate", "ribo_evidence", "n_aa",
                 "priority_score"]]
    return pd.concat([base, qc_table[qcols], dl, act[["actionability_score", "actionability_coverage"]],
                      mech[mcols]], axis=1)


def _row(feature, score, value, expectation, passed):
    return {"control": feature, "score": score, "value": value, "expectation": expectation,
            "pass": passed if passed is not None else np.nan}


def control_report(disease_loci, mech_by_locus_name):
    """Check the control expectations of LOGIC.md Step 5 that the current inputs can evaluate."""
    d = disease_loci
    g = lambda gene, col: d.loc[f"ctrl:{gene}", col] if f"ctrl:{gene}" in d.index else np.nan
    rows = [
        _row("HULC", "D1_HCC", f"log2FC {g('HULC', 'log2FC_paired'):.2f}, FDR {g('HULC', 'fdr_paired'):.2g}",
             "up in tumour, FDR < 0.05", g("HULC", "log2FC_paired") > 0 and g("HULC", "fdr_paired") < 0.05),
        _row("LINC00261", "D1_HCC",
             f"log2FC {g('LINC00261', 'log2FC_paired'):.2f}, FDR {g('LINC00261', 'fdr_paired'):.2g}",
             "down in tumour, FDR < 0.05", g("LINC00261", "log2FC_paired") < 0 and g("LINC00261", "fdr_paired") < 0.05),
        _row("COL1A1", "D1_chronic", f"rho_fibrosis {g('COL1A1', 'hoang_r'):.2f}, p {g('COL1A1', 'hoang_p'):.2g}, "
                                     f"D1_chronic {g('COL1A1', 'stat_D1_chronic'):.2f}",
             "positive association with fibrosis stage (rho > 0, p < 0.05)",
             g("COL1A1", "hoang_r") > 0 and g("COL1A1", "hoang_p") < 0.05),
    ]
    for gene in ("BIRC5", "CYTOR"):
        rows.append(_row(gene, "D2", f"HR_OS {g(gene, 'HR_OS'):.2f}, FDR {g(gene, 'fdr_OS'):.2g}", "HR > 1",
                         g(gene, "HR_OS") > 1))
    rows += [
        _row("PVT1", "D3", f"rho_CN {g('PVT1', 'cn_rho'):.2f}, FDR {g('PVT1', 'cn_fdr'):.2g}, gain "
                           f"{g('PVT1', 'cn_gain_freq'):.2f}, D3 {g('PVT1', 'stat_D3'):.2f}",
             "rho_CN > 0 (FDR < 0.05), D3 > 0",
             g("PVT1", "cn_rho") > 0 and g("PVT1", "cn_fdr") < 0.05 and g("PVT1", "stat_D3") > 0),
    ]
    for gene, note in (("MEG3", ""), ("GSTP1", "; chosen post hoc")):
        rows.append(_row(gene, "D3", f"rho_meth {g(gene, 'meth_rho'):.2f} (raw {g(gene, 'meth_rho_raw'):.2f}), "
                                     f"FDR {g(gene, 'meth_fdr'):.2g}",
                         f"rho_meth < 0 (FDR < 0.05){note}", g(gene, "meth_rho") < 0 and g(gene, "meth_fdr") < 0.05))
    for gene in ("PNPLA3", "TM6SF2"):
        rows.append(_row(gene, "D4", f"max -log10 p {g(gene, 'stat_D4'):.0f} ({g(gene, 'gwas_top_snp')})",
                         "liver-trait GWAS association", g(gene, "stat_D4") >= -np.log10(5e-8)))
    rows += [
        _row("HULC", "D5", f"percentile {g('HULC', 'pct_D5'):.2f}", "high", g("HULC", "pct_D5") >= 0.9),
        _row("LINC01090", "D5", f"liver papers {g('LINC01090', 'n_liver'):.0f}",
             "no LINC01089 papers credited", g("LINC01090", "n_liver") <= 1),
        _row("HNF1A-AS1", "D6", f"D6 {g('HNF1A-AS1', 'comp_D6'):.2f}", "liver-specific (D6 >= 0.75)",
             g("HNF1A-AS1", "comp_D6") >= 0.75),
    ]
    for col in ("disease_score", "disease_score_fitted"):
        if col not in d:
            continue
        bg90 = d.loc[~d["is_control"], col].quantile(0.9)
        for gene in ("HULC", "H19"):
            rows.append(_row(gene, col, f"{g(gene, col):.3f} (90th pct {bg90:.3f})",
                             "> 90th percentile of the background", g(gene, col) > bg90))
    rows.append(_row("GPC3", "actionability", f"tumour specificity {g('GPC3', 'tumour_specificity'):.2f}",
                     "tumour-specific (> 0); HLA binding not evaluated", g("GPC3", "tumour_specificity") > 0))

    expected = {"M1": ["MTLN", "SPAAR"], "M2": ["GAS5", "ZFAS1"], "M3": ["HNF1A-AS1", "LINC00261"],
                "M4": ["MALAT1", "NEAT1"], "M5": ["NORAD"], "M6": ["SNHG1", "MIR22HG"]}
    classes = [f"M{i}" for i in range(1, 7)]
    for cls, genes in expected.items():
        for gene in genes:
            if gene not in mech_by_locus_name.index:
                rows.append(_row(gene, cls, "not a lincorfz locus", "top class" if cls != "M2" else "high", None))
                continue
            sc = mech_by_locus_name.loc[gene, [f"{c}_score" for c in classes]].astype(float)
            top = classes[int(np.argmax(sc.to_numpy()))]
            val = f"{cls} {sc[f'{cls}_score']:.2f}; top {top}"
            ok = sc[f"{cls}_score"] >= 0.6 if cls == "M2" else top == cls
            rows.append(_row(gene, cls, val, "top class" if cls != "M2" else "high (>= 0.6)", ok))
    return pd.DataFrame(rows)
