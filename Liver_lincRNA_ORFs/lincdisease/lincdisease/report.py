"""
Final tables (LOGIC.md, "Outputs") and the control panel check (Step 5).

Per ORF: QC flags and tier, copy-attribution confidence, the locus' disease components and score, the
actionability score, mechanism class scores, standardized class scores and dominant class.
Actionability (per ORF) = equal-weight composite of the percentiles (among all ORFs, each weighted by 1 / number
of ORFs of its locus) of three parts: the locus' tumour specificity, the number of lincORF-specific HLA-I
immunopeptides, and the predicted HLA-I presentation (MHCflurry; distinct binder peptides over the HepG2, Huh7 and
common alleles; missing where the peptide step has not been run).
"""
import operator

import numpy as np
import pandas as pd

from . import controls, scoring

ACTIONABILITY = ["tumour_specificity", "n_hla_specific", "hla_binding"]


def actionability(orfs, disease_loci, peptides, hla_table=None):
    """Per ORF: the ACTIONABILITY parts (act_), their percentiles (act_pct_) and the composite score."""
    f = pd.DataFrame(index=orfs.index)
    f["tumour_specificity"] = orfs["locus_id"].map(disease_loci["tumour_specificity"]).to_numpy()
    f["n_hla_specific"] = orfs["orf_id"].map(peptides["n_hla_specific"]).fillna(0).to_numpy()
    f["hla_binding"] = (orfs["orf_id"].map(hla_table["hla_any_n_binders"]).astype(float).to_numpy()
                        if hla_table is not None else np.nan)
    pct = scoring.component_percentiles(f, f, ACTIONABILITY, weights=scoring.locus_weights(orfs["locus_id"]))
    res = scoring.composite(pct, dict.fromkeys(ACTIONABILITY, 1.0))
    return pd.concat([f.add_prefix("act_"), pct.add_prefix("act_pct_"),
                      res.rename(columns={"score": "actionability_score", "coverage": "actionability_coverage"})],
                     axis=1)


def orf_table(orfs, qc_table, disease_loci, mech, act):
    """lincdisease.orfs.tsv: ORF identity, QC, the locus' disease columns, actionability and mechanism scores."""
    keep_disease = ["disease_score_noD5", "disease_score", "disease_coverage", "disease_direction",
                    "hcc_score", "hcc_coverage", "hcc_direction", "chronic_score", "chronic_coverage",
                    "chronic_direction",
                    "comp_D1", "comp_D2", "comp_D3", "comp_D4", "comp_D5", "comp_D6",
                    "log2FC_paired", "fdr_paired", "auc_vs_gtex", "fdr_vs_gtex", "hcc_rep_r", "hcc_rep_fdr",
                    "hcc_rep_pairs", "hcc_replicated",
                    "stat_D1_chronic", "chronic_fdr", "hoang_r", "hoang_p", "govaere_r", "govaere_p", "suppli_r",
                    "suppli_p",
                    "hcv_r", "hcv_p",
                    "D1_chronic_cell", "chronic_cell_fdr", "parent_gene", "parent_rho", "parent_driven",
                    "HR_OS", "fdr_OS", "HR_PFI",
                    "fdr_PFI", "gwas_top_trait", "gwas_top_snp", "gwas_disease_max", "gwas_enzyme_max",
                    "gwas_hcc_max", "gwas_chronic_max", "n_liver", "n_liver_cancer", "n_liver_chronic", "n_review",
                    "n_retracted", "tumour_specificity"]
    keep_disease += [c for c in ("disease_score_fitted",) if c in disease_loci]
    dl = disease_loci[keep_disease].reindex(orfs["locus_id"])
    dl.index = orfs.index
    mcols = [c for c in mech.columns if c.endswith(("_score", "_z", "_coverage")) or c == "dominant_class"]
    qcols = ["confidence_tier", "copy_attribution", "segdup_fracMatch", "grch38_false_dup", "antisense_to_coding",
             "small_RNA_overlap", "multimapper_dominated", "coding_gene_extension", "tcga_quantifiable"]
    base = orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate", "ribo_evidence", "n_aa",
                 "priority_score"]]
    return pd.concat([base, qc_table[qcols], dl, act[["actionability_score", "actionability_coverage"]],
                      mech[mcols]], axis=1)


SCREEN_COLS = ["crispr_screen_value", "prensner_value", "prensner_resolved", "prensner_hepg2_lfc", "chen_value",
               "chen_K562_phenotype", "chen_K562_fdr", "chen_iPSC_phenotype", "chen_iPSC_fdr",
               "phylop_orf_minus_nonorf"]
TIER_ORDER = {"clean": 0, "caution": 1, "expressed_paralog": 2, "likely_artefact": 3}


def screen_supported(table, support):
    """ORFs with a viability phenotype in a published CRISPR screen (screen value >= 1), candidates or not."""
    s = support.set_index("orf_id")[SCREEN_COLS].reindex(table["orf_id"]).set_axis(table.index)
    out = pd.concat([table, s], axis=1)
    out = out[out["crispr_screen_value"] >= 1]
    key = out["confidence_tier"].map(TIER_ORDER)
    return (out.assign(_tier=key)
            .sort_values(["crispr_screen_value", "_tier", "disease_score_noD5"], ascending=[False, True, False])
            .drop(columns="_tier"))


OPERATORS = {">": operator.gt, ">=": operator.ge, "<": operator.lt, "<=": operator.le, "==": operator.eq}


def _format(v):
    return f"{v:.3g}" if isinstance(v, (float, np.floating)) else str(v)


def _missing(v):
    return v is None or (not isinstance(v, str) and pd.isna(v))


def control_report(disease_loci, mech_by_locus_name, hla_control=None):
    """
    Evaluate the checks of controls.PANEL (LOGIC.md Step 5): one row per check with the values used, the
    expectation (from the check's conditions) and pass (True / False; missing = not evaluated: the control or
    one of its values is missing).
    """
    d = disease_loci
    classes = [f"M{i}" for i in range(1, 7)]

    def values(check):
        """The check's source row as a dict, or a reason why there is none."""
        if check.source == "disease":
            key = f"ctrl:{check.gene}"
            return d.loc[key].to_dict() if key in d.index else "not quantified"
        if check.source == "mechanism":
            if check.gene not in mech_by_locus_name.index:
                return "not a lincorfz locus"
            row = mech_by_locus_name.loc[check.gene].to_dict()
            z = [float(row[f"{c}_z"]) for c in classes]
            row["top_class"] = classes[int(np.nanargmax(z))] if not np.all(np.isnan(z)) else np.nan
            return row
        return hla_control.iloc[0].to_dict() if hla_control is not None and len(hla_control) else "not run"

    rows = []
    for check in controls.PANEL:
        row = values(check)
        conds, texts = [], []
        for col, op, threshold in check.conditions:
            text = f"{col} {op} {_format(threshold)}"
            if threshold == controls.BACKGROUND_90:
                threshold = d.loc[~d["is_control"], col].quantile(0.9) if col in d else np.nan
                text = f"{col} {op} {_format(threshold)}, the {controls.BACKGROUND_90}"
            conds.append((col, op, threshold))
            texts.append(text)
        expectation = f"{check.expectation} ({'; '.join(texts)})"
        if isinstance(row, str):
            rows.append({"control": check.gene, "score": check.score, "value": row, "expectation": expectation,
                         "pass": np.nan})
            continue
        cols = list(dict.fromkeys([c for c, _, _ in conds] + list(check.show)))
        value = ", ".join(f"{c} {_format(row.get(c, np.nan))}" for c in cols)
        if any(_missing(row.get(c)) or _missing(t) for c, _, t in conds):
            passed = np.nan
        else:
            passed = all(OPERATORS[op](row[c], t) for c, op, t in conds)
        rows.append({"control": check.gene, "score": check.score, "value": value, "expectation": expectation,
                     "pass": passed})
    return pd.DataFrame(rows)
