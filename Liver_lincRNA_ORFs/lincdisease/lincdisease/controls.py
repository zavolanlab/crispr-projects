"""
Positive / negative controls (LOGIC.md, Step 5). They are GENCODE v42 genes processed by the same code as the
candidates: quantified as features 'ctrl:<gene>' and scored against the same background. The mechanism controls
check that the features of a class are defined as intended (e.g. a divergent lncRNA has the divergent-position
feature), not that the mechanism is detected: several of them pass on annotation features alone.

PANEL is the one definition of the checks: report.control_report evaluates its conditions and writes the
expectation text from them, so the two cannot disagree.
"""
from typing import NamedTuple

BACKGROUND_90 = "90th percentile of the background"   # threshold: the column's 90th percentile over the loci


class Check(NamedTuple):
    """
    One control check. Every condition (column, operator, threshold) must hold; a missing value of any of them
    makes the check 'not evaluated'. `show`: further columns printed with the value. `source`: the table the
    columns come from: disease (disease.loci.tsv, row ctrl:<gene>), mechanism (mechanism.orfs.tsv, the best-ranked
    ORF of the locus named <gene>; top_class = the class with the highest z) or hla (hla.control.tsv).
    """
    score: str
    gene: str
    expectation: str
    conditions: tuple
    show: tuple = ()
    source: str = "disease"


GWAS_SIGNIFICANT = 7.3      # -log10(5e-8)


def _top_class(cls, gene):
    return Check(cls, gene, "top class", (("top_class", "==", cls),), (f"{cls}_score", f"{cls}_z"), "mechanism")


PANEL = [
    Check("D1_HCC", "HULC", "up in tumour", (("log2FC_paired", ">", 0), ("fdr_paired", "<", 0.05))),
    Check("D1_HCC", "LINC00261", "down in tumour", (("log2FC_paired", "<", 0), ("fdr_paired", "<", 0.05))),
    Check("D1_HCC replication", "HULC", "up in tumour in the independent cohorts",
          (("hcc_rep_r", ">", 0), ("hcc_rep_fdr", "<", 0.05)), ("hcc_rep_pairs",)),
    Check("D1_HCC replication", "LINC00261", "down in tumour in the independent cohorts",
          (("hcc_rep_r", "<", 0), ("hcc_rep_fdr", "<", 0.05)), ("hcc_rep_pairs",)),
    Check("D1_chronic", "COL1A1", "positive association with fibrosis stage (Hoang)",
          (("hoang_r", ">", 0), ("hoang_p", "<", 0.05)), ("stat_D1_chronic",)),
    Check("D2", "BIRC5", "worse overall survival with higher expression", (("HR_OS", ">", 1),), ("fdr_OS",)),
    Check("D2", "CYTOR", "worse overall survival with higher expression", (("HR_OS", ">", 1),), ("fdr_OS",)),
    Check("D3", "PVT1", "expression follows copy number",
          (("cn_rho", ">", 0), ("cn_fdr", "<", 0.05), ("stat_D3", ">", 0)), ("cn_gain_freq",)),
    Check("D4", "PNPLA3", "genome-wide significant liver-trait GWAS association",
          (("stat_D4", ">=", GWAS_SIGNIFICANT),), ("gwas_top_snp",)),
    Check("D4", "TM6SF2", "genome-wide significant liver-trait GWAS association",
          (("stat_D4", ">=", GWAS_SIGNIFICANT),), ("gwas_top_snp",)),
    Check("D5", "HULC", "high", (("pct_D5", ">=", 0.9),)),
    Check("D5", "LINC01090", "no LINC01089 papers credited", (("n_liver", "<=", 1),)),
    Check("D6", "HNF1A-AS1", "liver-specific", (("comp_D6", ">=", 0.75),)),
    *[Check(col, gene, "high", ((col, ">", BACKGROUND_90),))
      for col, genes in (("disease_score", ("HULC", "H19")), ("disease_score_fitted", ("HULC", "H19")),
                         ("hcc_score", ("HULC", "H19")), ("chronic_score", ("COL1A1", "PNPLA3", "TM6SF2")))
      for gene in genes],
    Check("actionability", "GPC3", "tumour-specific", (("tumour_specificity", ">", 0),)),
    Check("actionability", "GPC3", "GPC3 144-152 predicted HLA-I binder", (("presentation_percentile", "<=", 2),),
          ("peptide", "allele", "affinity"), "hla"),
    _top_class("M1", "MTLN"), _top_class("M1", "SPAAR"),
    Check("M2", "GAS5", "high", (("M2_score", ">=", 0.6),), ("M2_z", "top_class"), "mechanism"),
    Check("M2", "ZFAS1", "high", (("M2_score", ">=", 0.6),), ("M2_z", "top_class"), "mechanism"),
    _top_class("M3", "HNF1A-AS1"), _top_class("M3", "LINC00261"),
    _top_class("M4", "MALAT1"), _top_class("M4", "NEAT1"),
    _top_class("M5", "NORAD"),
    _top_class("M6", "SNHG1"), _top_class("M6", "MIR22HG"),
]

# neighbours needed by the M3 controls (cis co-expression)
NEIGHBOURS = {"HNF1A-AS1": "HNF1A", "LINC00261": "FOXA2"}


def control_genes():
    genes = {c.gene for c in PANEL} | set(NEIGHBOURS.values())
    return sorted(genes)
