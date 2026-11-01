"""
D4 germline genetics (LOGIC.md, Step 3): GWAS Catalog associations with liver traits near each locus.

Liver traits are whitelisted by MAPPED_TRAIT (EFO labels): liver diseases (HCC, MASLD, fibrosis, cirrhosis,
alcoholic liver disease, viral hepatitis, ...) and the liver-enzyme traits ALT, AST and GGT. The LD window of
a locus is approximated by its genomic span +- `flank` (default 50 kb). Statistic D4 = the largest
-log10(p) of a whitelisted association in the window (0 without one), reported separately for disease and
enzyme traits.
"""
import zipfile

import numpy as np
import pandas as pd

LIVER_DISEASE_TRAITS = {
    "hepatocellular carcinoma", "metabolic dysfunction-associated steatotic liver disease",
    "non-alcoholic fatty liver disease", "nonalcoholic steatohepatitis", "liver fibrosis measurement",
    "cirrhosis of liver", "alcoholic liver disease", "liver disorder", "liver fat measurement",
    "chronic hepatitis B virus infection", "chronic hepatitis C virus infection", "hepatitis B virus infection",
    "hepatitis C virus infection", "susceptibility to hepatitis B infection measurement",
    "drug-induced liver injury", "Abnormality of the liver", "biliary liver cirrhosis",
}
HCC_TRAITS = {"hepatocellular carcinoma"}     # the rest of LIVER_DISEASE_TRAITS is chronic liver disease
LIVER_ENZYME_TRAITS = {
    "serum alanine aminotransferase amount", "alanine aminotransferase 1 measurement",
    "aspartate aminotransferase measurement", "serum gamma-glutamyl transferase measurement",
    "aspartate aminotransferase to alanine aminotransferase ratio",
}


def read_liver_associations(zip_path):
    """Whitelisted associations: chrom ('chr..'), pos (1-based), mlog10p, trait, trait_class, snp, pubmed."""
    with zipfile.ZipFile(zip_path) as z:
        with z.open(z.namelist()[0]) as fh:
            df = pd.read_csv(fh, sep="\t", dtype=str, usecols=["CHR_ID", "CHR_POS", "PVALUE_MLOG", "MAPPED_TRAIT",
                                                                  "SNPS", "PUBMEDID"])
    df = df.dropna(subset=["CHR_ID", "CHR_POS", "MAPPED_TRAIT"])
    traits = df["MAPPED_TRAIT"].str.split(", ")
    df = df.assign(trait=traits).explode("trait")
    df["trait_class"] = np.where(df["trait"].isin(LIVER_DISEASE_TRAITS), "disease",
                                 np.where(df["trait"].isin(LIVER_ENZYME_TRAITS), "enzyme", ""))
    df = df[df["trait_class"] != ""]
    # multi-SNP rows ('1;1', '123;456') are rare; keep the first position
    pos = pd.to_numeric(df["CHR_POS"].str.split(r"[;x ]").str[0], errors="coerce")
    out = pd.DataFrame({"chrom": "chr" + df["CHR_ID"].str.split(r"[;x ]").str[0], "pos": pos,
                        "mlog10p": pd.to_numeric(df["PVALUE_MLOG"], errors="coerce"),
                        "trait": df["trait"], "trait_class": df["trait_class"], "snp": df["SNPS"],
                        "pubmed": df["PUBMEDID"]})
    return out.dropna(subset=["pos", "mlog10p"]).astype({"pos": int})


def d4(spans, assoc, flank=50_000):
    """`spans`: DataFrame index feature_id with chrom/start0/end1. Max -log10 p in span +- flank per trait class."""
    rows = []
    by_chrom = {c: g.sort_values("pos") for c, g in assoc.groupby("chrom")}
    for fid, r in spans.iterrows():
        g = by_chrom.get(r["chrom"])
        rec = {"feature_id": fid, "gwas_disease_max": 0.0, "gwas_hcc_max": 0.0, "gwas_chronic_max": 0.0,
               "gwas_enzyme_max": 0.0, "gwas_n": 0,
               "gwas_top_trait": "", "gwas_top_snp": ""}
        if g is not None:
            lo, hi = r["start0"] - flank, r["end1"] + flank
            pos = g["pos"].to_numpy()
            w = g.iloc[np.searchsorted(pos, lo, "left"):np.searchsorted(pos, hi, "right")]
            if len(w):
                rec["gwas_n"] = len(w)
                for cls in ("disease", "enzyme"):
                    m = w.loc[w["trait_class"] == cls, "mlog10p"]
                    rec[f"gwas_{cls}_max"] = float(m.max()) if len(m) else 0.0
                hcc = w["trait"].isin(HCC_TRAITS)
                for name, sel in (("hcc", hcc), ("chronic", (w["trait_class"] == "disease") & ~hcc)):
                    m = w.loc[sel, "mlog10p"]
                    rec[f"gwas_{name}_max"] = float(m.max()) if len(m) else 0.0
                top = w.sort_values("mlog10p").iloc[-1]
                rec["gwas_top_trait"], rec["gwas_top_snp"] = top["trait"], top["snp"]
        rows.append(rec)
    out = pd.DataFrame(rows).set_index("feature_id")
    out["D4"] = out[["gwas_disease_max", "gwas_enzyme_max"]].max(axis=1)
    # per-score parts: HCC; chronic liver disease (non-HCC liver disease traits and liver enzymes)
    out["D4_hcc"] = out["gwas_hcc_max"]
    out["D4_chronic"] = out[["gwas_chronic_max", "gwas_enzyme_max"]].max(axis=1)
    return out
