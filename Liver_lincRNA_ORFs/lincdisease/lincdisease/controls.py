"""
Positive / negative controls (LOGIC.md, Step 5). They are GENCODE v42 genes processed by the same code as the
candidates: quantified as features 'ctrl:<gene>' and scored against the same background.
"""

# score -> list of (gene, expectation)
PANEL = {
    "D1_HCC": [("HULC", "up in tumour, FDR < 0.05"), ("LINC00261", "down in tumour, FDR < 0.05")],
    "D1_chronic": [("COL1A1", "positive association with fibrosis stage")],
    "D2": [("BIRC5", "HR > 1"), ("CYTOR", "HR > 1")],
    "D3": [("PVT1", "positive CN-expression correlation"), ("MEG3", "methylation-expression anticorrelation")],
    "D4": [("PNPLA3", "liver-trait GWAS association"), ("TM6SF2", "liver-trait GWAS association")],
    "D5": [("HULC", "high"), ("LINC01090", "no LINC01089 papers credited")],
    "D6": [("HNF1A-AS1", "liver-specific")],
    "disease_score": [("HULC", "> 90th percentile"), ("H19", "> 90th percentile")],
    "actionability": [("GPC3", "GPC3 144-152 predicted binder, tumour-specific")],
    "M1": [("MTLN", "top class"), ("SPAAR", "top class")],
    "M2": [("GAS5", "high"), ("ZFAS1", "high")],
    "M3": [("HNF1A-AS1", "top class"), ("LINC00261", "top class")],
    "M4": [("MALAT1", "top class"), ("NEAT1", "top class")],
    "M5": [("NORAD", "top class")],
    "M6": [("SNHG1", "top class"), ("MIR22HG", "top class")],
}

# neighbours needed by the M3 controls (cis co-expression)
NEIGHBOURS = {"HNF1A-AS1": "HNF1A", "LINC00261": "FOXA2"}


def control_genes():
    genes = {g for entries in PANEL.values() for g, _ in entries} | set(NEIGHBOURS.values())
    return sorted(genes)
