"""
M1 (encoded peptide): hits of published CRISPR knock-out screens of non-canonical ORFs (LOGIC.md, Step 4).

Chen et al. 2020 (Science 367:1140, "Pervasive functional translation of non-canonical human open reading frames";
Tables S1, S3, S4): Cas9 screens of ~2,300 ORFs in K562 and iPSCs (read_chen, chen_hits). Their published hit
calls are not in the tables and their stated procedure (non-targeting null) makes nearly every ORF a hit, as Cas9
cutting alone costs ~0.015 in phenotype; hits are therefore re-called against the cutting-shifted null.

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
CHEN_CUTTING = -0.015       # mean phenotype from Cas9 cutting alone (Chen et al. 2020, Methods)


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
    out["hepg2_lfc"] = out["name"].map(lfc)
    out["secondary"] = out["name"].map(s32.drop_duplicates("Gene").set_index("Gene")["Phenotype type"])
    out["value"] = screen_value(out["screened"], out["phenotype"], out["secondary"])
    out["resolved"] = (out["phenotype"].astype(bool) & out["secondary"].notna()).astype(float)
    return out


def screen_value(screened, phenotype, secondary):
    """2 ORF-specific phenotype, 1 unresolved phenotype, 0 no or non-ORF-specific phenotype, NaN not screened."""
    v = np.where(phenotype.astype(bool), np.where(secondary.isna(), 1.0, secondary.isin(ORF_SPECIFIC) * 2.0), 0.0)
    return pd.Series(np.where(screened, v, np.nan), index=screened.index)


def match(orfs, screen, prefix, report=None):
    """Per ORF: matched screen ORF names, the best (maximum) screen value and each `report` column ({column: "min" or
    "max"}) aggregated over them (NaN without a match)."""
    names = peptide.same_stop(orfs, screen)
    s = screen.set_index("name")
    best = lambda col, f: names.map(lambda n: f(s.loc[n.split(";"), col]) if n else np.nan)
    out = {f"{prefix}_orfs": names, f"{prefix}_value": best("value", lambda x: x.max())}
    for c, how in (report or {}).items():
        out[f"{prefix}_{c}"] = best(c, lambda x: getattr(x.astype(float), how)())
    return pd.DataFrame(out, index=orfs.index)


# ---------------------------------------------------------------- hg19 -> hg38 (UCSC chain file, points only)

def read_chain(path):
    """UCSC chain file -> per target chrom: arrays of ungapped blocks (t_start, t_end, q_chrom, q_start, q_strand,
    q_size, score)."""
    import gzip
    blocks = {}
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for ln in fh:
            f = ln.split()
            if not f:
                continue
            if f[0] == "chain":
                score, t_chrom, t_pos = float(f[1]), f[2], int(f[5])
                q_chrom, q_size, q_strand, q_pos = f[7], int(f[8]), f[9], int(f[10])
                rows = blocks.setdefault(t_chrom, [])
                continue
            size = int(f[0])
            rows.append((t_pos, t_pos + size, q_chrom, q_pos, q_strand, q_size, score))
            if len(f) == 3:
                t_pos += size + int(f[1])
                q_pos += size + int(f[2])
    return {c: pd.DataFrame(r, columns=["t_start", "t_end", "q_chrom", "q_start", "q_strand", "q_size", "score"])
            for c, r in blocks.items()}


def lift_points(chain, chrom, pos0, strand):
    """0-based positions -> (chrom, pos0, strand) in the query assembly via the best-scoring chain block containing
    the position; (None, nan, None) if unmapped."""
    out = []
    for c, p, s in zip(chrom, pos0, strand):
        b = chain.get(c)
        hit = None if b is None else b[(b["t_start"] <= p) & (b["t_end"] > p)]
        if hit is None or hit.empty:
            out.append((None, np.nan, None))
            continue
        h = hit.loc[hit["score"].idxmax()]
        q = h["q_start"] + (p - h["t_start"])
        if h["q_strand"] == "-":
            q, s = h["q_size"] - 1 - q, {"+": "-", "-": "+"}[s]
        out.append((h["q_chrom"], int(q), s))
    return pd.DataFrame(out, columns=["chrom", "pos0", "strand"])


# ---------------------------------------------------------------- Chen et al. 2020

def top3_phenotype(gammas):
    """ORF phenotype of Chen et al.: mean of the 3 sgRNA phenotypes with the largest magnitude."""
    g = np.asarray(gammas, dtype=float)
    return float(g[np.argsort(-np.abs(g))[:3]].mean())


def chen_hits(guides, neg_gamma, cutting=CHEN_CUTTING, fdr=0.05, min_effect=-0.02):
    """
    Per ORF (guides: orf/gamma): top-3 phenotype, one-sided Mann-Whitney p of its sgRNAs against the non-targeting
    sgRNAs shifted by the mean Cas9 cutting phenotype, BH-FDR, hit = FDR < fdr and phenotype <= min_effect.
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    null = np.asarray(neg_gamma, dtype=float) + cutting
    rows = [(o, top3_phenotype(g), mannwhitneyu(g, null, alternative="less").pvalue)
            for o, g in guides.groupby("orf")["gamma"]]
    r = pd.DataFrame(rows, columns=["name", "phenotype", "p"]).set_index("name")
    r["fdr"] = multipletests(r["p"], method="fdr_bh")[1]
    r["hit"] = (r["fdr"] < fdr) & (r["phenotype"] <= min_effect)
    return r


def read_chen(d, chain):
    """
    Chen et al. 2020 tables in directory d (aay0262_table_s1/s3/s4.xlsx): ORFs with their hg38 stop codon (lifted
    from hg19; S1 coordinates are 0-based, stop = one past the stop codon in the direction of transcription; blocks
    hold only the stop codon, for same_stop), K562 / iPSC phenotypes and hits (chen_hits), screen value (1 = hit in
    either cell line, 0 = screened without a hit).
    """
    t = lambda n: f"{d}/aay0262_table_{n}.xlsx"
    s1 = pd.concat(pd.read_excel(t("s1"), sheet_name=None, header=3).values()).drop_duplicates("ORF name")
    s3 = pd.read_excel(t("s3"), header=2)
    s4 = pd.read_excel(t("s4"), sheet_name=None, header=None, skiprows=5)
    neg = s3["gene"] == "negative_control"
    out = None
    for cell in ("K562", "iPSC"):
        gamma = s4[f"{cell} sgRNA phenotype"].set_index(0)[3]
        g = pd.DataFrame({"orf": s3["ORF name"], "gamma": s3["guide name"].map(gamma), "neg": neg}).dropna()
        r = chen_hits(g[~g["neg"]], g.loc[g["neg"], "gamma"]).add_prefix(f"{cell}_")
        out = r if out is None else out.join(r, how="outer")
    c = s1.set_index("ORF name").reindex(out.index)
    last = c["stop genomic coordinate"] + np.where(c["strand"] == "+", -1, 1)
    lifted = lift_points(chain, c["chromosome"], last, c["strand"])
    out = out.reset_index().assign(chrom=lifted["chrom"].to_numpy(), strand=lifted["strand"].to_numpy(),
                                   orf_type=c["ORF type"].to_numpy())
    end = lifted["pos0"].to_numpy() + np.where(out["strand"] == "+", 1, 0)
    out["blocks"] = [[(int(e) - 3, int(e))] if s == "+" else [(int(e), int(e) + 3)] if s == "-" else []
                     for e, s in zip(end, out["strand"])]
    out = out[out["chrom"].notna()].reset_index(drop=True)
    out["value"] = (out["K562_hit"].fillna(False).astype(bool) | out["iPSC_hit"].fillna(False).astype(bool)).astype(float)
    return out


def combined_value(prensner_value, prensner_resolved, chen_value):
    """CRISPR screen feature: Prensner's secondary-screen verdict where it exists (the only test of ORF
    specificity), otherwise the larger of the two screen values (missing if neither screened the ORF)."""
    v = pd.concat([prensner_value, chen_value], axis=1).max(axis=1)
    return v.where(~prensner_resolved, prensner_value)
