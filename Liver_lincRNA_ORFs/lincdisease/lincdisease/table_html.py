"""
Self-contained HTML table of the candidate and screen-supported ORFs (`lincdisease.table.html`): sortable columns
ordered from the summary scores to the detail, column groups that can be hidden, a filter by set, and an index
of the summary scores, D1-D6 and M1-M6. No external resources.
"""
import html
import json

import numpy as np
import pandas as pd

# (group, column, label, kind, description); kind: text, score (0-1, shaded), num, p (p-value / FDR), int, bool
COLUMNS = [
    ("ORF", "locus_name", "Locus", "text", "lincRNA locus (gene symbol, Ensembl or StringTie id)"),
    ("ORF", "orf_label", "ORF", "text", "ORF label from lincorfz"),
    ("ORF", "set", "Set", "text", "candidate (curated list of 46 ORFs), screen (CRISPR screen value >= 1), or both"),
    ("ORF", "confidence_tier", "QC tier", "text",
     "clean / caution / expressed_paralog / likely_artefact (Methods 6)"),
    ("ORF", "n_aa", "Length (aa)", "int", "ORF length in amino acids"),
    ("ORF", "ribo_evidence", "Ribo-seq", "text", "translation evidence class from lincorfz"),

    ("Summary scores", "disease_score_noD5", "Disease score w/o D5", "score",
     "disease score without prior literature (D5): the headline score, unbiased by how well studied a locus is"),
    ("Summary scores", "hcc_score", "HCC score", "score",
     "D1 HCC 0.25, D2 0.25, D3 0.20, D4 HCC 0.10, D5 HCC 0.10, D6 0.10; 0.5 = typical lincRNA"),
    ("Summary scores", "chronic_score", "Chronic score", "score",
     "D1 chronic 0.50, D4 chronic 0.20, D5 chronic 0.15, D6 0.15; 0.5 = typical lincRNA"),
    ("Summary scores", "disease_score", "Disease score", "score", "D1 0.25, D2 0.25, D3 0.20, D4-D6 0.10 each"),
    ("Summary scores", "actionability_score", "Actionability", "score",
     "tumour specificity, predicted HLA-I presentation, HLA-I immunopeptides (1/3 each); per ORF"),
    ("Summary scores", "crispr_screen_value", "CRISPR screen", "num",
     "2 = ORF-specific viability phenotype (Prensner 2021), 1 = phenotype (unresolved or Chen 2020 hit), 0 = none"),
    ("Summary scores", "dominant_class", "Dominant class", "text",
     "mechanism class with the highest standardized score z, above the median (z > 0) and >= 0.5 ahead of the next, "
     "else none"),
    ("Summary scores", "hcc_direction", "HCC direction", "text", "from D1 HCC, D2, D3 (informative: percentile >= 0.75)"),
    ("Summary scores", "chronic_direction", "Chronic direction", "text", "up / down with chronic liver disease (D1 chronic)"),

    ("Disease components", "comp_D1", "D1", "score", "disease-associated expression: 0.7 HCC + 0.3 chronic (percentile)"),
    ("Disease components", "comp_D2", "D2", "score", "prognosis, TCGA-LIHC Cox OS and PFI (percentile of |z|)"),
    ("Disease components", "comp_D3", "D3", "score", "somatic copy number (percentile)"),
    ("Disease components", "comp_D4", "D4", "score", "liver-trait GWAS within +-50 kb (percentile)"),
    ("Disease components", "comp_D5", "D5", "score", "liver literature, Europe PMC (percentile)"),
    ("Disease components", "comp_D6", "D6", "score", "liver context: GTEx liver level and specificity, HepG2/Huh7 (percentile)"),
    ("Disease components", "disease_direction", "Disease direction", "text", "oncogene-like / suppressor-like / inconsistent / none"),
    ("Disease components", "hcc_coverage", "HCC coverage", "score", "fraction of the HCC score weight with data"),
    ("Disease components", "chronic_coverage", "Chronic coverage", "score", "fraction of the chronic score weight with data"),
    ("Disease components", "disease_coverage", "Disease coverage", "score", "fraction of the disease score weight with data"),

    *[("Mechanism classes", f"M{i}_score", f"M{i}", "score", d) for i, d in enumerate([
        "encoded peptide", "translation-coupled RNA regulation", "cis-regulation of a neighbour",
        "trans-acting RNA, nuclear", "trans-acting RNA, cytoplasmic", "small-RNA host"], 1)],
    *[("Mechanism z", f"M{i}_z", f"M{i} z", "num",
       f"M{i} standardized score: normal quantile of the class-score percentile among all ORFs (locus-weighted)")
      for i in range(1, 7)],

    ("D1 HCC detail", "log2FC_paired", "log2FC T/N", "num", "mean log2 fold change, tumour vs adjacent normal (50 pairs)"),
    ("D1 HCC detail", "fdr_paired", "FDR T/N", "p", "Wilcoxon signed-rank, BH"),
    ("D1 HCC detail", "auc_vs_gtex", "AUC vs GTEx", "num", "tumours vs GTEx liver, (r + 1) / 2 (reported, not scored: study batch)"),
    ("D1 HCC detail", "fdr_vs_gtex", "FDR vs GTEx", "p", "Mann-Whitney, BH"),
    ("D1 HCC detail", "hcc_rep_r", "r replication", "num",
     "paired rank-biserial r, tumour vs adjacent normal in 4 independent cohorts (GSE77314, GSE77509, GSE94660, "
     "GSE56545), pair-weighted mean"),
    ("D1 HCC detail", "hcc_rep_fdr", "FDR replication", "p", "Stouffer over the 4 cohorts, BH"),
    ("D1 HCC detail", "hcc_replicated", "Replicated", "bool", "same direction as D1 HCC and replication FDR < 0.05"),

    ("D1 chronic detail", "stat_D1_chronic", "D1 chronic", "num",
     "weighted mean of cohort partial rank correlations (fibrosis 2/3: Hoang + Govaere by samples tested; HCV 1/3)"),
    ("D1 chronic detail", "chronic_fdr", "Chronic FDR", "p", "Stouffer, BH"),
    ("D1 chronic detail", "hoang_r", "r fibrosis", "num", "Hoang 2019, partial rank correlation with fibrosis stage"),
    ("D1 chronic detail", "govaere_r", "r fibrosis (Govaere)", "num",
     "Govaere 2020, partial rank correlation with fibrosis stage"),
    ("D1 chronic detail", "hoang_p", "p fibrosis", "p", ""),
    ("D1 chronic detail", "hcv_r", "r HCV", "num", "Boldanova 2017, HCV vs non-HCV"),
    ("D1 chronic detail", "hcv_p", "p HCV", "p", ""),
    ("D1 chronic detail", "suppli_r", "r NAFL/NASH", "num", "Suppli 2019 (reported, not scored: technical confounding)"),
    ("D1 chronic detail", "suppli_p", "p NAFL/NASH", "p", ""),
    ("D1 chronic detail", "D1_chronic_cell", "D1 chronic, cell-adjusted", "num",
     "additionally adjusted for liver cell-type marker scores (reported, not scored)"),
    ("D1 chronic detail", "chronic_cell_fdr", "FDR cell-adjusted", "p", ""),
    ("D1 chronic detail", "parent_gene", "Parent gene", "text", "homology parent (Ensembl id)"),
    ("D1 chronic detail", "parent_rho", "rho parent", "num", "median Spearman correlation with the parent over cohorts"),
    ("D1 chronic detail", "parent_driven", "Parent-driven", "bool",
     "signal follows a disease-associated parent; statistic re-computed with the parent as covariate"),

    ("D2 detail", "HR_OS", "HR OS", "num", "hazard ratio per SD of expression, overall survival (adjusted for age, stage, log size factor)"),
    ("D2 detail", "fdr_OS", "FDR OS", "p", ""),
    ("D2 detail", "HR_PFI", "HR PFI", "num", "progression-free interval"),
    ("D2 detail", "fdr_PFI", "FDR PFI", "p", ""),

    ("D3 detail", "stat_D3", "D3", "num", "signed somatic statistic (> 0: expression follows a copy gain)"),
    ("D3 detail", "cn_rho", "rho CN", "num", "purity-adjusted Spearman, copy number vs expression"),
    ("D3 detail", "cn_fdr", "FDR CN", "p", ""),
    ("D3 detail", "cn_gain_freq", "CN gain freq", "num", "fraction of tumours with log2 ratio > 0.3"),
    ("D3 detail", "cn_loss_freq", "CN loss freq", "num", "fraction of tumours with log2 ratio < -0.3"),

    ("D4 detail", "gwas_top_trait", "Top GWAS trait", "text", "liver trait of the strongest association within +-50 kb"),
    ("D4 detail", "gwas_top_snp", "Top SNP", "text", ""),
    ("D4 detail", "gwas_hcc_max", "-log10 p HCC", "num", "hepatocellular carcinoma"),
    ("D4 detail", "gwas_chronic_max", "-log10 p chronic", "num", "chronic liver disease traits and liver enzymes"),
    ("D4 detail", "gwas_enzyme_max", "-log10 p enzymes", "num", "ALT, AST, GGT"),

    ("D5 detail", "n_liver", "Liver papers", "int", "Europe PMC: symbol AND RNA term AND liver term"),
    ("D5 detail", "n_liver_cancer", "Cancer papers", "int", "liver papers with a cancer term"),
    ("D5 detail", "n_liver_chronic", "Chronic papers", "int", "liver papers with a chronic liver disease term"),
    ("D5 detail", "n_review", "Reviews", "int", ""),
    ("D5 detail", "n_retracted", "Retracted", "int", ""),

    ("Actionability detail", "tumour_specificity", "Tumour specificity", "num",
     "median log2 ncov in TCGA-LIHC tumours minus the largest GTEx tissue median (testis excluded)"),
    ("Actionability detail", "hla_any_n_binders", "HLA-I binders", "int",
     "distinct 8-11-mers with MHCflurry presentation percentile <= 2 (HepG2, Huh7, common alleles)"),
    ("Actionability detail", "n_hla_specific", "HLA-I peptides", "int", "lincORF-specific immunopeptides identified"),
    ("Actionability detail", "n_specific_peptides", "Proteomics peptides", "int",
     "proteomics peptides absent from canonical proteins"),

    ("CRISPR screen detail", "prensner_value", "Prensner", "num", "Prensner 2021: 2 ORF-specific, 1 unresolved, 0 none"),
    ("CRISPR screen detail", "prensner_hepg2_lfc", "Prensner HepG2 LFC", "num", "median sgRNA log2 fold change in HepG2"),
    ("CRISPR screen detail", "chen_value", "Chen", "num", "Chen 2020: 1 hit (K562 or iPSC, re-called), 0 none"),
    ("CRISPR screen detail", "chen_K562_phenotype", "Chen K562", "num", "phenotype (gamma) in K562"),
    ("CRISPR screen detail", "chen_K562_fdr", "Chen K562 FDR", "p", ""),
    ("CRISPR screen detail", "chen_iPSC_phenotype", "Chen iPSC", "num", "phenotype (gamma) in iPSC"),
    ("CRISPR screen detail", "chen_iPSC_fdr", "Chen iPSC FDR", "p", ""),
    ("CRISPR screen detail", "phylop_orf_minus_nonorf", "PhyloP ORF - exon", "num",
     "mean PhyloP100way over the ORF minus over the rest of the transcript's exons"),

    *[("Mechanism features (percentiles)", f"pct_{f}", lab, "score", desc) for f, lab, desc in [
        ("m1_translation_weighted", "M1 translation",
         "pooled log10 Bayes factor of translation (within length bin, weighted by copy attribution)"),
        ("m1_orf_share_weighted", "M1 ORF P-site share",
         "fraction of the transcript's P-sites in the ORF (within length bin, weighted by copy attribution)"),
        ("m1_peptides", "M1 peptides", "identified peptides (proteomics + HLA-I)"),
        ("m1_independent_translation", "M1 GENCODE Ribo-seq", "same stop as a GENCODE Ribo-seq ORF (+ peptide support)"),
        ("m1_phylocsf_weighted", "M1 PhyloCSF",
         "within length and HepG2-expression bin, weighted by copy-attribution confidence"),
        ("m1_phylop_contrast_weighted", "M1 PhyloP contrast", "ORF minus exon PhyloP, weighted by copy attribution"),
        ("m1_sequence_features_weighted", "M1 sequence",
         "TM / signal-like / Pfam, and order (within length bin, weighted by copy attribution)"),
        ("m1_crispr_screen", "M1 CRISPR", "screen value"),
        ("m2_nmd", "M2 NMD", "stop > 50 nt upstream of the last exon junction"),
        ("m2_nmd_sensitivity", "M2 NMD sensitivity",
         "mean log2FC upon UPF1 and SMG6 + SMG7 knock-down (HeLa; within HepG2-expression bin)"),
        ("m3_neighbour_position", "M3 position", "divergent or antisense to a protein-coding gene"),
        ("m3_coexpression", "M3 co-expression", "excess co-expression with the nearest coding gene"),
        ("m3_enhancer_chromatin", "M3 enhancer", "enhancer-like chromatin at the TSS"),
        ("m3_low_polysome", "M3 low polysome", "low mean ribosome load (within HepG2-expression bin)"),
        ("m4_abundance", "M4 abundance", "HepG2 log2 TPM"),
        ("m4_eclip_nuclear", "M4 eCLIP nuclear", "binding by nuclear RBPs"),
        ("m5_cytoplasmic", "M5 cytoplasmic", "HepG2 cytoplasmic enrichment"),
        ("m5_eclip_cytoplasmic", "M5 eCLIP cytoplasmic", "binding by cytoplasmic RBPs"),
        ("m6_small_rna_host", "M6 small RNAs", "snoRNA / scaRNA / miRNA genes hosted"),
    ]],

    ("QC detail", "copy_attribution", "Copy attribution", "score",
     "confidence that the evidence belongs to this copy (1 outside duplications)"),
    ("QC detail", "segdup_fracMatch", "Segdup identity", "num", "identity of the overlapping segmental duplication"),
    ("QC detail", "grch38_false_dup", "GRCh38 false dup", "bool", ""),
    ("QC detail", "antisense_to_coding", "Antisense to coding", "bool", ""),
    ("QC detail", "small_RNA_overlap", "Small-RNA overlap", "bool", ""),
    ("QC detail", "multimapper_dominated", "Multimapper-dominated", "bool", ""),
    ("QC detail", "coding_gene_extension", "Coding-gene extension", "bool", ""),
    ("QC detail", "tcga_quantifiable", "Quantifiable", "bool", ">= 100 quantifiable positions in unstranded RNA-seq"),

    ("Identifiers", "orf_id", "ORF id", "text", ""),
    ("Identifiers", "locus_id", "Locus id", "text", ""),
    ("Identifiers", "priority_score", "lincorfz priority", "num", "lincorfz ranking score"),
]

INDEX = [
    ("Summary scores", [
        ("HCC score", "Relevance to hepatocellular carcinoma: weighted composite of D1 HCC, D2, D3, D4 HCC, D5 HCC "
                      "and D6 percentiles."),
        ("Chronic score", "Relevance to chronic liver disease (MASLD/MASH, fibrosis, viral hepatitis): D1 chronic, "
                          "D4 chronic, D5 chronic and D6."),
        ("Disease score", "Overall liver disease relevance, mostly HCC (chronic disease expression carries 7.5%)."),
        ("Actionability", "Per ORF: tumour-specific expression, predicted HLA-I presentation and identified HLA-I "
                          "peptides; suitability as a tumour antigen."),
        ("Reading the scores", "Every score is a percentile composite against the 2,057 expressed lincRNA loci (or "
                               "11,910 ORFs) of the run: 0.5 is a typical lincRNA, 1 the most extreme. Missing "
                               "components count 0.5, so a score near 0.5 with low coverage means unknown."),
    ]),
    ("Disease components (per locus)", [
        ("D1", "Disease-associated expression. HCC: rank-biserial effect of TCGA-LIHC tumours vs paired adjacent "
               "normals. Chronic: partial rank correlation with fibrosis stage (Hoang and Govaere) and HCV status "
               "(Boldanova), adjusted for technical covariates and parent genes."),
        ("D2", "Prognosis: Cox models of overall survival and progression-free interval in TCGA-LIHC, adjusted for "
               "age, stage and tumour composition (log size factor)."),
        ("D3", "Somatic copy number: frequent copy gains or losses that expression follows in TCGA-LIHC, "
               "purity-adjusted."),
        ("D4", "Germline genetics: strongest liver-trait GWAS association within 50 kb."),
        ("D5", "Prior knowledge: liver papers in Europe PMC (exact symbol, RNA term, liver term)."),
        ("D6", "Liver context: expression level and specificity in GTEx liver, and in HepG2/Huh7."),
    ]),
    ("Mechanism classes (per ORF)", [
        ("M1", "Encoded peptide: translation strength, peptides, conservation, sequence features, CRISPR "
               "screens; features a paralogous copy can supply are weighted by copy attribution, length-dependent "
               "ones ranked within ORF-length bins."),
        ("M2", "Translation-coupled RNA regulation: NMD substrate, measured NMD sensitivity."),
        ("M3", "Cis-regulation of a neighbour: divergent/antisense position, excess co-expression with the nearest "
               "coding gene, enhancer-like chromatin, low polysome association."),
        ("M4", "Trans-acting nuclear RNA: abundance, binding by nuclear RBPs."),
        ("M5", "Trans-acting cytoplasmic RNA: cytoplasmic enrichment, binding by cytoplasmic RBPs."),
        ("M6", "Small-RNA host: snoRNA, scaRNA or miRNA genes in the locus."),
        ("Scores and dominance", "Class score = mean of its feature percentiles (absolute evidence; classes are "
                                 "not exclusive). z = normal quantile of the class score's percentile among all ORFs "
                                 "(each locus counts once), so classes with few or binary features do not win by "
                                 "their wider spread. Dominant class: highest z, above 0 and >= 0.5 ahead of the "
                                 "next."),
    ]),
]


def rows(orfs, support, mechanism, disease, peptides, hla, min_screen=1):
    """Candidates and screen-supported ORFs (screen value >= min_screen) with all COLUMNS."""
    t = orfs.set_index("orf_id", drop=False)
    t = t.join(support.set_index("orf_id"), how="left", suffixes=("", "_support"))
    t = t.join(mechanism.set_index("orf_id")[[c for c in mechanism.columns if c.startswith("pct_") and c not in t]])
    t = t.join(peptides[["n_specific_peptides", "n_hla_specific"]]).join(hla[["hla_any_n_binders"]])
    locus_cols = [c for c in ("stat_D3", "cn_rho", "cn_fdr", "cn_gain_freq", "cn_loss_freq") if c in disease]
    t = t.join(disease[locus_cols], on="locus_id").copy()          # copy: the joins leave a fragmented frame
    cand = t["is_candidate"].astype(bool)
    screen = t["crispr_screen_value"].ge(min_screen)
    t["set"] = np.select([cand & screen, cand], ["both", "candidate"], "screen")
    return t[cand | screen].sort_values("disease_score_noD5", ascending=False)


def _value(v, kind):
    """A table value as JSON for the page: None when missing, else bool / float / str by column kind."""
    if not isinstance(v, str) and pd.isna(v):
        return None
    if kind == "bool":
        return bool(v) if not isinstance(v, str) else v == "True"
    if kind in ("score", "num", "p", "int"):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None
    return str(v)


def build(table, title="lincdisease: candidate and screen-supported ORFs", note=""):
    """The HTML page: COLUMNS present in `table` (rows()), the INDEX and the summary line filled into TEMPLATE."""
    cols = [c for c in COLUMNS if c[1] in table.columns]
    data = [[_value(r[c[1]], c[3]) for c in cols] for _, r in table.iterrows()]
    meta = [{"group": g, "col": c, "label": lab, "kind": k, "desc": d} for g, c, lab, k, d in cols]
    index_html = "".join(
        f"<section><h3>{html.escape(h)}</h3><dl>"
        + "".join(f"<dt>{html.escape(a)}</dt><dd>{html.escape(b)}</dd>" for a, b in items)
        + "</dl></section>" for h, items in INDEX)
    counts = table["set"].value_counts().to_dict()
    summary = (f"{len(table)} ORFs: {counts.get('candidate', 0) + counts.get('both', 0)} candidates, "
               f"{counts.get('screen', 0) + counts.get('both', 0)} with a CRISPR screen phenotype "
               f"({counts.get('both', 0)} both).")
    return (TEMPLATE.replace("__TITLE__", html.escape(title))
            .replace("__SUMMARY__", html.escape(summary + (" " + note if note else "")))
            .replace("__INDEX__", index_html)
            .replace("__META__", json.dumps(meta))
            .replace("__DATA__", json.dumps(data, allow_nan=False)))


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>lincdisease ORF table</title>
<style>
:root {
  --bg: #ffffff; --fg: #1d2430; --muted: #5d6878; --line: #d9dee5; --head: #f2f5f8; --head2: #e6ebf1;
  --accent: #2e5a88; --hover: #f7f9fb; --chip: #eef2f6; --heat: 46, 110, 170;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #14181e; --fg: #e3e8ee; --muted: #9aa5b3; --line: #2c333d; --head: #1c222a; --head2: #232a33;
    --accent: #8db7e2; --hover: #1a2027; --chip: #232a33; --heat: 96, 160, 225;
  }
}
:root[data-theme="dark"] {
  --bg: #14181e; --fg: #e3e8ee; --muted: #9aa5b3; --line: #2c333d; --head: #1c222a; --head2: #232a33;
  --accent: #8db7e2; --hover: #1a2027; --chip: #232a33; --heat: 96, 160, 225;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--fg);
       font: 14px/1.45 -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }
header, .index, .controls, footer { padding: 0 16px; max-width: 1400px; }
h1 { font-size: 20px; margin: 18px 0 4px; }
.summary { color: var(--muted); margin: 0 0 12px; }
details.index-box { margin: 0 16px 12px; border: 1px solid var(--line); border-radius: 6px; max-width: 1368px; }
details.index-box > summary { cursor: pointer; padding: 8px 12px; font-weight: 600; }
.index { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 4px 24px; padding-bottom: 8px; }
.index h3 { font-size: 14px; margin: 8px 0 4px; color: var(--accent); }
.index dl { margin: 0; }
.index dt { font-weight: 600; margin-top: 6px; }
.index dd { margin: 0 0 0 0; color: var(--muted); }
.controls { display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: center; margin-bottom: 10px; }
.controls input[type=search] { padding: 5px 8px; border: 1px solid var(--line); border-radius: 4px;
                               background: var(--bg); color: var(--fg); min-width: 200px; }
.seg button { border: 1px solid var(--line); background: var(--bg); color: var(--fg); padding: 4px 10px;
              cursor: pointer; font: inherit; }
.seg button:first-child { border-radius: 4px 0 0 4px; } .seg button:last-child { border-radius: 0 4px 4px 0; }
.seg button + button { border-left: none; }
.seg button.on { background: var(--accent); color: var(--bg); border-color: var(--accent); }
details.groups { position: relative; }
details.groups > summary { cursor: pointer; color: var(--accent); }
.group-list { position: absolute; z-index: 5; background: var(--bg); border: 1px solid var(--line);
              border-radius: 6px; padding: 8px 12px; margin-top: 4px; box-shadow: 0 4px 16px rgba(0,0,0,.12);
              column-count: 2; column-gap: 24px; width: max-content; max-width: 90vw; }
.group-list label { display: block; white-space: nowrap; }
.wrap { overflow: auto; max-height: calc(100vh - 120px); border-top: 1px solid var(--line);
        border-bottom: 1px solid var(--line); }
table { border-collapse: separate; border-spacing: 0; font-size: 13px; }
th, td { padding: 4px 8px; border-bottom: 1px solid var(--line); border-right: 1px solid var(--line);
         white-space: nowrap; background: var(--bg); }
thead th { position: sticky; z-index: 2; background: var(--head); }
thead tr.g th { top: 0; background: var(--head2); text-align: left; font-weight: 600; color: var(--accent); }
thead tr.c th { top: 27px; cursor: pointer; user-select: none; font-weight: 600; }
thead tr.c th:hover { color: var(--accent); }
thead tr.c th .arrow { color: var(--muted); font-size: 11px; }
th.first, td.first { position: sticky; left: 0; z-index: 1; font-weight: 600; }
thead th.first { z-index: 3; }
tbody tr:hover td { background: var(--hover); }
td.n { text-align: right; font-variant-numeric: tabular-nums; }
td.na { color: var(--muted); text-align: center; }
td.score { text-align: right; font-variant-numeric: tabular-nums; }
.chip { background: var(--chip); border-radius: 10px; padding: 1px 7px; font-size: 12px; }
footer { color: var(--muted); font-size: 12px; margin: 10px 0 20px; }
@media (max-width: 600px) { .index { grid-template-columns: 1fr; } .group-list { column-count: 1; } }
</style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <p class="summary">__SUMMARY__</p>
</header>
<details class="index-box"><summary>Index: summary scores, D1-D6, M1-M6</summary><div class="index">__INDEX__</div></details>
<div class="controls">
  <div class="seg" id="sets">
    <button data-set="all" class="on">All</button><button data-set="candidate">Candidates</button><button data-set="screen">Screen-supported</button>
  </div>
  <input type="search" id="q" placeholder="Search locus, ORF, trait...">
  <details class="groups"><summary>Column groups</summary><div class="group-list" id="groups"></div></details>
  <span class="summary" id="shown"></span>
</div>
<div class="wrap"><table id="t"><thead></thead><tbody></tbody></table></div>
<footer>Click a column header to sort (again to reverse); hover a header for its definition. Score columns are
percentiles against the background (0.5 = typical lincRNA) and are shaded by value. Missing values: &ndash;.</footer>
<script>
const META = __META__;
const DATA = __DATA__;
const groups = [...new Set(META.map(m => m.group))];
const hidden = new Set(["Mechanism z", "Mechanism features (percentiles)", "QC detail", "Identifiers"]);
let sortCol = META.findIndex(m => m.col === "disease_score_noD5"), sortDir = -1, setFilter = "all", query = "";
const setIdx = META.findIndex(m => m.col === "set");

function fmt(v, kind) {
  if (v === null) return "–";
  if (kind === "bool") return v ? "yes" : "no";
  if (kind === "int") return Math.round(v).toString();
  if (kind === "score") return v.toFixed(2);
  if (kind === "p") return v === 0 ? "0" : (v < 0.001 ? v.toExponential(1) : v.toPrecision(2));
  if (kind === "num") { const a = Math.abs(v); return a !== 0 && (a < 0.01 || a >= 1e4) ? v.toExponential(2) : (+v.toPrecision(3)).toString(); }
  return v;
}
function visibleCols() { return META.map((m, i) => i).filter(i => !hidden.has(META[i].group)); }
function render() {
  const cols = visibleCols();
  const thead = document.querySelector("#t thead");
  let g = '<tr class="g">', last = null, span = 0;
  const spans = [];
  cols.forEach(i => { if (META[i].group !== last) { if (last !== null) spans.push([last, span]); last = META[i].group; span = 0; } span++; });
  if (last !== null) spans.push([last, span]);
  spans.forEach(([name, n], k) => { g += `<th colspan="${n}"${k === 0 ? ' class="first"' : ''}>${name}</th>`; });
  g += "</tr>";
  let c = '<tr class="c">';
  cols.forEach((i, k) => {
    const m = META[i];
    const arrow = i === sortCol ? (sortDir < 0 ? " ▼" : " ▲") : "";
    c += `<th data-i="${i}"${k === 0 ? ' class="first"' : ''} title="${(m.col + (m.desc ? ": " + m.desc : "")).replace(/"/g, "&quot;")}">${m.label}<span class="arrow">${arrow}</span></th>`;
  });
  thead.innerHTML = g + c + "</tr>";
  thead.querySelectorAll("tr.c th").forEach(th => th.onclick = () => {
    const i = +th.dataset.i;
    if (i === sortCol) sortDir = -sortDir; else { sortCol = i; sortDir = META[i].kind === "text" ? 1 : -1; }
    render();
  });
  const q = query.toLowerCase();
  let rows = DATA.filter(r => (setFilter === "all" || r[setIdx] === setFilter || r[setIdx] === "both") &&
    (!q || r.some(v => v !== null && String(v).toLowerCase().includes(q))));
  rows.sort((a, b) => {
    const x = a[sortCol], y = b[sortCol];
    if (x === null && y === null) return 0; if (x === null) return 1; if (y === null) return -1;
    return (x < y ? -1 : x > y ? 1 : 0) * sortDir;
  });
  const heat = getComputedStyle(document.documentElement).getPropertyValue("--heat").trim();
  document.querySelector("#t tbody").innerHTML = rows.map(r => "<tr>" + cols.map((i, k) => {
    const m = META[i], v = r[i], first = k === 0 ? " first" : "";
    if (v === null) return `<td class="na${first}">–</td>`;
    if (m.kind === "score") return `<td class="score${first}" style="background: rgba(${heat}, ${(Math.max(0, Math.min(1, v)) * 0.45).toFixed(3)})">${fmt(v, m.kind)}</td>`;
    if (m.col === "set") return `<td class="${first}"><span class="chip">${v}</span></td>`;
    const cls = m.kind === "text" ? "" : "n";
    return `<td class="${cls}${first}">${fmt(v, m.kind)}</td>`;
  }).join("") + "</tr>").join("");
  document.getElementById("shown").textContent = `${rows.length} of ${DATA.length} ORFs shown`;
}
document.getElementById("groups").innerHTML = groups.map(g =>
  `<label><input type="checkbox" data-g="${g}"${hidden.has(g) ? "" : " checked"}${g === "ORF" ? " disabled" : ""}> ${g}</label>`).join("");
document.querySelectorAll("#groups input").forEach(cb => cb.onchange = () => {
  cb.checked ? hidden.delete(cb.dataset.g) : hidden.add(cb.dataset.g); render(); });
document.querySelectorAll("#sets button").forEach(b => b.onclick = () => {
  document.querySelectorAll("#sets button").forEach(x => x.classList.remove("on")); b.classList.add("on");
  setFilter = b.dataset.set; render(); });
document.getElementById("q").oninput = e => { query = e.target.value; render(); };
render();
</script>
</body>
</html>
"""
