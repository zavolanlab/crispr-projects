"""
M1 (encoded peptide): hits of published CRISPR knock-out screens of non-canonical ORFs (LOGIC.md, Step 4).

Prensner et al. 2021 (Nat Biotechnol 39:697, "Noncanonical open reading frames encode functional proteins essential
for cancer cell survival"; supplementary tables, NIHMS1698134):
* Table S2: the ORFeome (553 ORFs, GRCh38 exon coordinates, `CRISPR phenotype` = viability phenotype in the primary
  Cas9 screen of 8 cell lines, among them HepG2);
* Table S22: primary screen, log2 fold change of every sgRNA (day 21 vs plasmid) per cell line;
* Table S32: secondary tiling screen of the hits, classifying whether the phenotype is specific to the ORF.

A lincORF matches an ORFeome ORF with the same stop codon (as for the GENCODE Ribo-seq catalog, peptide.same_stop).
Screen value: 2 = viability phenotype specific to the ORF (secondary screen: "Specific to ORF" or "Specific to ORF
and transcript subregions"); 1 = viability phenotype in the primary screen, not resolved by the secondary screen;
0 = screened without phenotype, or a phenotype that the secondary screen attributes to the transcript, the genomic
region or an adjacent coding gene (Cas9 cutting there does not implicate the peptide); missing = not screened.
"""
import numpy as np
import pandas as pd

from . import peptide

ORF_SPECIFIC = {"Specific to ORF", "Specific to ORF and transcript subregions"}


def _coords(s):
    return [int(x) for x in str(s).replace(" ", "").split(",") if x]


def read_prensner(tables_1_17, tables_21_32):
    """ORFeome ORFs with blocks (start0, end1), primary-screen phenotype, HepG2 LFC and secondary-screen class."""
    s2 = pd.read_excel(tables_1_17, sheet_name="Table S2", header=4)
    s22 = pd.read_excel(tables_21_32, sheet_name="Table S22", header=2)
    s32 = pd.read_excel(tables_21_32, sheet_name="Table S32", header=3)
    rows = []
    for r in s2.itertuples(index=False):
        a, b = _coords(r[s2.columns.get_loc("Exon Start (GCh38)")]), _coords(r[s2.columns.get_loc("Exon End (GCh38)")])
        blocks = sorted((min(x, y) - 1, max(x, y)) for x, y in zip(a, b))
        rows.append((r.Chr, r.Strand, r.Name, blocks))
    out = pd.DataFrame(rows, columns=["chrom", "strand", "name", "blocks"])
    out["phenotype"] = s2["CRISPR phenotype"].fillna(0).astype(int).to_numpy()
    lfc = s22.groupby("Gene Symbol")["LFC HepG2_pDNA_D21"].median()
    out["screened"] = out["name"].isin(lfc.index)
    out["hepg2_median_lfc"] = out["name"].map(lfc)
    out["secondary"] = out["name"].map(s32.drop_duplicates("Gene").set_index("Gene")["Phenotype type"])
    out["value"] = screen_value(out["screened"], out["phenotype"], out["secondary"])
    return out


def screen_value(screened, phenotype, secondary):
    """2 ORF-specific phenotype, 1 unresolved phenotype, 0 no or non-ORF-specific phenotype, NaN not screened."""
    v = np.where(phenotype.astype(bool), np.where(secondary.isna(), 1.0, secondary.isin(ORF_SPECIFIC) * 2.0), 0.0)
    return pd.Series(np.where(screened, v, np.nan), index=screened.index)


def match(orfs, screen, prefix):
    """Per ORF: matched screen ORF names and the best screen value / HepG2 LFC among them (NaN without a match)."""
    names = peptide.same_stop(orfs, screen)
    s = screen.set_index("name")
    best = lambda col, f: names.map(lambda n: f(s.loc[n.split(";"), col]) if n else np.nan)
    return pd.DataFrame({f"{prefix}_orfs": names,
                         f"{prefix}_value": best("value", lambda x: x.max()),
                         f"{prefix}_hepg2_lfc": best("hepg2_median_lfc", lambda x: x.min())}, index=orfs.index)
