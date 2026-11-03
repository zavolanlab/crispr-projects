"""
D4 germline genetics, GTEx v8 liver eQTLs (LOGIC.md, Step 3), for features with a GENCODE gene: lincRNA loci
through their GENCODE ids (name_source_gene_ids, lncRNA_gene_ids), controls through their gene name. GTEx v8
uses GENCODE v26; ids are matched without version.

Per feature (over its genes): eGene at q < 0.05 (GTEx `qval`), lead variant, its -log10 p and slope, number of
significant variant-gene pairs, and the link to a liver-trait GWAS association (genetics.py whitelist):
  eqtl_gwas_shared  a significant eQTL variant is a GWAS hit (same chromosome and position)
  eqtl_gwas_near    the lead eQTL variant lies within `proxy_distance` (10 kb) of a GWAS hit (LD proxy; no LD
                    reference panel is used)
  eqtl_gwas_mlog10p largest -log10 p of a GWAS hit satisfying either link (0 if none)
Features without a GENCODE gene, or whose genes GTEx does not test, are missing.
"""
import numpy as np
import pandas as pd

GTEX_TAR = "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar"
EGENES = "GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz"
PAIRS = "GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz"
QVAL = 0.05
PROXY_DISTANCE = 10_000


def _ids(value):
    if not isinstance(value, str):
        return []
    return [v.strip().split(".")[0] for v in value.split(",") if v.strip().startswith("ENSG")]


def feature_genes(loci, control_names, egenes):
    """DataFrame feature_id/gene_id (unversioned) linking loci and controls to GTEx genes."""
    rows = [(fid, g) for fid, r in loci.iterrows()
            for g in sorted(set(_ids(r.get("name_source_gene_ids")) + _ids(r.get("lncRNA_gene_ids"))))]
    by_name = egenes.assign(g=egenes["gene_id"].str.split(".").str[0]).groupby("gene_name")["g"].first()
    rows += [(f"ctrl:{n}", by_name[n]) for n in control_names if n in by_name.index]
    return pd.DataFrame(rows, columns=["feature_id", "gene_id"]).drop_duplicates()


def parse_variant(v):
    """GTEx variant ids chr_pos_ref_alt_b38 -> (chrom, pos)."""
    p = v.str.split("_")
    return p.str[0], p.str[1].astype(int)


def per_feature(links, egenes, pairs, gwas, proxy_distance=PROXY_DISTANCE, qval=QVAL):
    """
    The eQTL columns of the module docstring per feature. `links`: feature_genes; `egenes`, `pairs`: GTEx liver
    eGene and significant variant-gene pair tables; `gwas`: genetics.read_liver_associations.
    """
    eg = egenes.assign(g=egenes["gene_id"].str.split(".").str[0]).set_index("g")
    pr = pairs.assign(g=pairs["gene_id"].str.split(".").str[0])
    pr["chrom"], pr["pos"] = parse_variant(pr["variant_id"])
    gw = gwas.groupby(["chrom", "pos"])["mlog10p"].max()
    gw_by_chrom = {c: g.sort_values("pos") for c, g in gwas.groupby("chrom")}
    rows = []
    for fid, gl in links.groupby("feature_id"):
        genes = [g for g in gl["gene_id"] if g in eg.index]
        if not genes:
            continue
        e = eg.loc[genes]
        lead = e.sort_values("pval_nominal").iloc[0]
        sig = e[e["qval"] < qval]
        p = pr[pr["g"].isin(sig.index)]
        rec = {"feature_id": fid, "eqtl_genes": ",".join(genes), "eqtl_egene": len(sig) > 0,
               "eqtl_lead_gene": lead.name, "eqtl_lead_variant": lead["variant_id"],
               "eqtl_lead_rsid": lead.get("rs_id_dbSNP151_GRCh38p7", ""),
               "eqtl_lead_mlog10p": -np.log10(lead["pval_nominal"]), "eqtl_lead_slope": lead["slope"],
               "eqtl_qval": lead["qval"], "eqtl_n_signif_pairs": len(p),
               "eqtl_gwas_shared": False, "eqtl_gwas_near": False, "eqtl_gwas_mlog10p": 0.0}
        if len(sig):
            best = 0.0
            shared = [gw[k] for k in zip(p["chrom"], p["pos"]) if k in gw.index]
            if shared:
                rec["eqtl_gwas_shared"], best = True, max(shared)
            lsig = sig.sort_values("pval_nominal").iloc[0]
            c, pos = lsig["variant_id"].split("_")[:2]
            g = gw_by_chrom.get(c)
            if g is not None:
                w = g[(g["pos"] - int(pos)).abs() <= proxy_distance]
                if len(w):
                    rec["eqtl_gwas_near"], best = True, max(best, w["mlog10p"].max())
            rec["eqtl_gwas_mlog10p"] = best
        rows.append(rec)
    return pd.DataFrame(rows).set_index("feature_id")
