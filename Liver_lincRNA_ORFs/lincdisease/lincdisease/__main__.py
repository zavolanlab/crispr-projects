"""
Command-line driver of lincdisease; one sub-command per step of LOGIC.md, each writing into --out-dir. The workflow
(workflow/Snakefile) runs them; to run one by hand, from the output directory:

  PYTHONPATH=<lincdisease dir> python -m lincdisease --out-dir . <step> [options]

--out-dir comes before the step. Options of the steps that start from the lincorfz run:
  --results-dir R (lincorfz results/ directory, the background), --candidates T (candidate ORF table),
  --data-dir D (download cache of public resources).

  fetch      --data-dir D                       all downloads (login node)
  qc, homology, peptide, support, regions, literature, genetics, somatic, encode, mechanism, disease, weights,
  report     R T D [step options]
  manifest   --data-dir D                       recount3 samples of TCGA-LIHC and GTEx (login node)
  quant      [--sample-set quant] [--worker i --n-workers n] [--features ... --sums-dir sums_<name>]
                                                stream recount3 bigWigs (login node)
  scaling    --data-dir D                       size factors of the TCGA-LIHC / GTEx samples
  chronic    --manifest | --data-dir D [--govaere-dir ...]     D1 chronic liver disease (quant_sra/)
  nmd        --manifest | --data-dir D          M2 NMD sensitivity (quant_nmd/)
  hcc        --manifest | --data-dir D          D1 HCC replication (quant_hcc/)
  html                                          lincdisease.table.html

`python -m lincdisease --out-dir . <step> --help` lists the options of a step.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

from . import (chromatin, chronic, controls, disease, eclip, encode, eqtl, expression, genetics, hcc, hla, homology,
               inputs, integrate, literature, mechanism, nmd, peptide, qc, quant, report, resources, scoring, screens,
               somatic, table_html, weights)


def _common(p):
    p.add_argument("--results-dir", required=True, help="lincorfz results/ directory (background)")
    p.add_argument("--candidates", required=True, help="candidate ORF table (lincorfz format)")
    p.add_argument("--data-dir", required=True, help="download cache for public resources")


def read_peptides(pdir):
    """Specific and copy-specific peptides per ORF from the lincorfz proteomics and immunopeptidomics results."""
    read = lambda *parts: pd.read_csv(os.path.join(pdir, *parts), sep="\t")
    return mechanism.specific_peptides(
        read("proteomics", "analysis", "lincORF_peptides.tsv"), read("proteomics", "fasta", "lincORF_entries.tsv"),
        read("peptidomics", "analysis", "immunopeptides.tsv"), read("peptidomics", "fasta", "lincORF_entries.tsv"))


def cmd_qc(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    peptides = read_peptides(args.peptides_dir)
    peptides.to_csv(os.path.join(args.out_dir, "peptides.orfs.tsv"), sep="\t", index_label="orf_id")
    kw = dict(segdups=resources.read_segdups(args.data_dir), false_dups=resources.read_false_dups(args.data_dir),
              copy_peptides=peptides["copy_specific"])
    feats_path = os.path.join(args.out_dir, "quant", "features.tsv")
    if os.path.exists(feats_path):
        feats = pd.read_csv(feats_path, sep="\t", index_col=0)
        kw["quantifiable"] = orfs.set_index("orf_id")["locus_id"].map(feats["quantifiable"])
    hpath = os.path.join(args.out_dir, "homology.orfs.tsv")
    if os.path.exists(hpath):
        hom = pd.read_csv(hpath, sep="\t")
        kw["homology"] = hom.loc[hom["parent"].notna(), ["orf_id", "parent", "pident", "qcov"]]
    res = qc.run_qc(orfs, **kw)
    if os.path.exists(hpath):
        h = hom.set_index("orf_id").reindex(orfs["orf_id"])
        res["homolog_parent_coverage"] = h["parent_coverage"].to_numpy()
    out = pd.concat([orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate"]],
                     res.drop(columns="orf_id")], axis=1)
    path = os.path.join(args.out_dir, "qc.orfs.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"wrote {path}: {len(out)} ORFs")
    print("tiers, all ORFs:", out["confidence_tier"].value_counts().to_dict())
    print("tiers, candidates:", out.loc[out["is_candidate"], "confidence_tier"].value_counts().to_dict())


def cmd_homology(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    gencode = quant.read_gtf_exons(resources.fetch("gencode_gtf", args.data_dir))
    res = homology.run(orfs, gencode, args.data_dir, os.path.join(args.out_dir, "homology"), diamond=args.diamond,
                       threads=args.threads)
    path = os.path.join(args.out_dir, "homology.orfs.tsv")
    res.to_csv(path, sep="\t")
    print(f"wrote {path}: with homolog {res['parent'].notna().sum()}")


def cmd_peptide(args):
    """M1 sequence features (TM, signal-like, Pfam, disorder), GENCODE Ribo-seq ORF matches, MHCflurry HLA-I
    presentation. hmmsearch, metapredict and MHCflurry run in a separate environment (--tool-env)."""
    orfs = inputs.load(args.results_dir, args.candidates)
    env = os.path.join(args.tool_env, "bin")
    out = orfs[["orf_id"]].join(peptide.sequence_features(orfs))
    rd = os.path.join(args.data_dir, "riboseq_orfs")
    p2 = peptide.read_bed12(os.path.join(rd, "phase2", "Ribo-seq_ORFs.comprehensive.bed"))
    p1 = peptide.read_bed12(os.path.join(rd, "phase1", "Ribo-seq_ORFs.peptideins.bed"))
    out["gencode_riboseq_orfs"] = peptide.same_stop(orfs, p2)
    out["gencode_riboseq_orf"] = out["gencode_riboseq_orfs"] != ""
    out["gencode_riboseq_peptide_support"] = peptide.same_stop(orfs, p1) != ""
    seqs = orfs.set_index("orf_id")["aa_seq"].fillna("").str.rstrip("*")
    pf = peptide.pfam(seqs, os.path.join(env, "hmmsearch"), os.path.join(args.data_dir, "Pfam-A.hmm"), cpu=args.threads)
    out = out.join(pf, on="orf_id")
    out["pfam_n_domains"] = out["pfam_n_domains"].fillna(0).astype(int)
    out["pfam_domains"] = out["pfam_domains"].fillna("")
    out = out.join(peptide.disorder(seqs, os.path.join(env, "python")).add_prefix("pep_"), on="orf_id")
    out.to_csv(os.path.join(args.out_dir, "peptide_features.orfs.tsv"), sep="\t", index=False)
    print(f"peptide features: Ribo-seq catalog {out['gencode_riboseq_orf'].sum()}, Pfam {(out['pfam_n_domains'] > 0).sum()}")
    if args.skip_hla:
        return
    menv = {"MHCFLURRY_DATA_DIR": os.path.join(args.data_dir, "mhcflurry")}
    hd = os.path.join(args.data_dir, "homology")
    proteins = [os.path.join(hd, f) for f in ("gencode.v42.pc_translations.fa.gz", "UP000005640_9606.fasta.gz",
                                              "UP000005640_9606_additional.fasta.gz")]
    res = hla.run(orfs, os.path.join(env, "python"), os.path.join(args.out_dir, "hla.orfs.tsv"), env=menv,
                  protein_fastas=proteins)
    print(f"HLA: ORFs with >= 1 non-self binder {(res['hla_any_n_binders'] > 0).sum()} of {len(res)}; "
          f"self-peptide binders left out: {res['hla_n_self_binders'].sum()}")


def cmd_support(args):
    """M1 conservation contrast (PhyloP100way, ORF vs the rest of its transcript's exons) and published CRISPR
    screen hits (Prensner et al. 2021 supplementary tables in <data_dir>/crispr_screens)."""
    orfs = inputs.load(args.results_dir, args.candidates)
    tx_exons = quant.read_gtf_exons(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.gtf"))
    regions = pd.read_csv(os.path.join(args.out_dir, "quant", "regions.tsv"), sep="\t")
    cons = peptide.conservation_contrast(orfs, tx_exons, bigwig=args.phylop, threads=args.threads, regions=regions)
    d = os.path.join(args.data_dir, "crispr_screens")
    p = os.path.join(d, "NIHMS1698134-supplement-Supplementary_Tables_")
    prensner = screens.match(orfs, screens.read_prensner(p + "1-17.xlsx", p + "21-32.xlsx"), "prensner",
                             {"hepg2_lfc": "min", "resolved": "max"})
    chain = screens.read_chain(resources.fetch("hg19_to_hg38_chain", args.data_dir))
    chen = screens.match(orfs, screens.read_chen(d, chain), "chen",
                         {"K562_phenotype": "min", "K562_fdr": "min", "iPSC_phenotype": "min", "iPSC_fdr": "min"})
    out = pd.concat([orfs[["orf_id"]], cons, prensner, chen], axis=1)
    out["crispr_screen_value"] = screens.combined_value(out["prensner_value"], out["prensner_resolved"] == 1,
                                                        out["chen_value"])
    out.to_csv(os.path.join(args.out_dir, "support.orfs.tsv"), sep="\t", index=False)
    print(f"support: PhyloP contrast for {out['phylop_orf_minus_nonorf'].notna().sum()} ORFs; screen matches "
          f"Prensner {(out['prensner_orfs'] != '').sum()}, Chen {(out['chen_orfs'] != '').sum()}; combined values "
          f"{out['crispr_screen_value'].value_counts().to_dict()}")


def cmd_regions(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    loci = inputs.loci(orfs)
    linc = quant.read_gtf_exons(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.gtf"))
    gencode = quant.read_gtf_exons(resources.fetch("gencode_gtf", args.data_dir))
    regions, feats = quant.build_regions(linc, gencode, loci, controls.control_genes(),
                                         min_quant_bases=args.min_quant_bases)
    d = os.path.join(args.out_dir, "quant")
    os.makedirs(d, exist_ok=True)
    regions.to_csv(os.path.join(d, "regions.tsv"), sep="\t", index=False)
    feats.to_csv(os.path.join(d, "features.tsv"), sep="\t", index_label="feature_id")
    print(f"{len(feats)} features, {len(regions)} regions; not quantifiable: {(~feats['quantifiable']).sum()}")


TCGA_COLUMNS = ["tcga_barcode", "cgc_sample_sample_type", "gdc_cases.case_id", "gdc_cases.submitter_id"]
GTEX_COLUMNS = ["SUBJID", "SMTS", "SMTSD", "SEX", "AGE"]


def tcga_gtex_projects(data_dir):
    """(source, project) of the quant/ sample set: TCGA-LIHC, then every GTEx tissue."""
    projects = pd.read_csv(resources.fetch("gtex_projects", data_dir), sep="\t")["project"].unique()
    return [("tcga", "LIHC")] + [("gtex", p) for p in sorted(set(projects) - {"STUDY_NA"})]   # STUDY_NA: unknown


def cmd_manifest(args):
    """TCGA-LIHC (all samples), GTEx liver (all) and up to --gtex-per-tissue samples of every other GTEx tissue."""
    parts = [quant.manifest("tcga", "LIHC", TCGA_COLUMNS)]
    for _, project in tcga_gtex_projects(args.data_dir)[1:]:
        m = quant.manifest("gtex", project, GTEX_COLUMNS)
        if project != "LIVER" and len(m) > args.gtex_per_tissue:
            m = m.sample(args.gtex_per_tissue, random_state=args.seed)
        parts.append(m)
    man = pd.concat(parts, ignore_index=True)
    d = os.path.join(args.out_dir, "quant")
    os.makedirs(d, exist_ok=True)
    man.to_csv(os.path.join(d, "samples.tsv"), sep="\t", index=False)
    print(man.groupby(["source", "project"]).size().to_string())


def cmd_scaling(args):
    """Size factors of the TCGA-LIHC and GTEx samples (quant/size_factors.tsv), over all of them together, since
    tumours, adjacent normals and the GTEx tissues are compared with each other (D1, D6, tumour specificity)."""
    ncov, samples, _ = expression.load(args.out_dir)
    pc = expression.protein_coding_ids(resources.fetch("gencode_gtf", args.data_dir))
    projects = {(s, p) for s, p in zip(samples["source"], samples["project"])}
    genes = pd.concat([quant.read_gene_sums(resources.fetch(f"{s}_{p}_genes", args.data_dir,
                                                            url=quant.gene_sums_url(s, p)), s, samples)
                       for s, p in tcga_gtex_projects(args.data_dir) if (s, p) in projects], axis=1)
    missing = samples["sample_id"][~samples["sample_id"].isin(genes.columns)]
    if len(missing):
        raise ValueError(f"no gene sums for {len(missing)} samples ({list(missing[:3])})")
    sf = expression.size_factors(genes[genes.index.isin(pc)])
    path = os.path.join(args.out_dir, "quant", "size_factors.tsv")
    sf.rename("size_factor").to_csv(path, sep="\t", index_label="sample_id")
    print(f"wrote {path}; median by group:")
    grp = samples.set_index("sample_id")
    grp = np.where(grp["source"] == "tcga", "TCGA " + grp["cgc_sample_sample_type"].astype(str), "GTEx " + grp["project"])
    print(sf.groupby(pd.Series(grp, index=samples["sample_id"]).reindex(sf.index)).median().round(2).to_string())


def tcga_gtex_expression(out_dir, pseudocount):
    """Size-factor scaled ncov of the TCGA-LIHC / GTEx samples, log2(ncov + pseudocount), samples, features and
    the size factors."""
    ncov, samples, feats = expression.load(out_dir)
    sf = expression.read_size_factors(out_dir)
    ncov = expression.scaled(ncov, sf)
    samples = samples[samples["sample_id"].isin(ncov.columns)]
    return ncov, disease.log_expr(ncov, pseudocount), samples, feats, sf


def cmd_quant(args):
    regions = pd.read_csv(os.path.join(args.out_dir, "quant", "regions.tsv"), sep="\t")
    d = os.path.join(args.out_dir, args.sample_set)
    samples = pd.read_csv(os.path.join(d, "samples.tsv"), sep="\t", low_memory=False)
    if args.project:
        samples = samples[samples["project"].isin(args.project)]
    if args.features:
        regions = regions[regions["feature_id"].isin(args.features)]
    samples = samples.iloc[args.worker::args.n_workers]
    path = quant.run(samples, regions, os.path.join(d, args.sums_dir), worker=f"{args.worker}", scratch=args.scratch)
    print(f"done: {path}")


def feature_names(results_dir, candidates):
    """feature_id -> symbol for every locus (its locus_name) and control ('ctrl:<gene>')."""
    loci = inputs.loci(inputs.load(results_dir, candidates))
    names = loci["locus_name"].copy()
    for g in controls.control_genes():
        names[f"ctrl:{g}"] = g
    return names


def cmd_chronic(args):
    """D1 chronic liver disease from the SRA cohorts: manifest (needs network) or statistics on the quantified sums."""
    d = os.path.join(args.out_dir, "quant_sra")
    if args.manifest:
        os.makedirs(d, exist_ok=True)
        chronic.manifest().to_csv(os.path.join(d, "samples.tsv"), sep="\t", index=False)
        return
    ncov, samples, feats = expression.load(args.out_dir, "quant_sra")
    # loci that cannot be measured in unstranded RNA-seq are not tested (and stay out of the BH family)
    ncov = ncov.loc[ncov.index.isin(expression.quantifiable(feats))]
    genes = pd.concat([quant.read_gene_sums(resources.fetch(f"sra_{p}_genes", args.data_dir,
                                                            url=quant.gene_sums_url("sra", p)), "sra", samples)
                       for p in chronic.COHORTS.values()], axis=1)
    qc = pd.concat([pd.read_csv(resources.fetch(f"sra_{p}_qc", args.data_dir,
                                                url=quant._metadata_url("sra", p, "recount_qc")),
                                sep="\t", low_memory=False) for p in chronic.COHORTS.values()])
    qc.index = "sra:" + qc["external_id"].astype(str)
    gtf = quant.read_gtf_exons(resources.fetch("gencode_gtf", args.data_dir)).drop_duplicates("gene_id")
    gtf["gene_id"] = gtf["gene_id"].str.split(".").str[0]
    symbol_to_id = dict(zip(gtf["gene_name"], gtf["gene_id"]))
    # parent gene of each locus: the homology parent of its ORF with the highest identity
    orf_qc = pd.read_csv(os.path.join(args.out_dir, "qc.orfs.tsv"), sep="\t").dropna(subset=["homolog_parent"])
    best = orf_qc.sort_values("homolog_pident", ascending=False).drop_duplicates("locus_id").set_index("locus_id")
    parents = best["homolog_parent"].map(symbol_to_id).dropna()
    pheno = samples.set_index("sample_id")[["cohort", "fibrosis_stage", "nas", "group"]]
    if args.govaere_dir:
        regions = pd.read_csv(os.path.join(args.out_dir, "quant", "regions.tsv"), sep="\t")
        g_ncov, g_genes, g_qc, g_pheno = chronic.govaere(args.govaere_dir, args.govaere_runs, args.govaere_soft,
                                                         regions, exclude=args.govaere_exclude)
        print(f"Govaere: {len(g_pheno)} runs")
        missing = ncov.index.difference(g_ncov.index)
        if len(missing):
            raise ValueError(f"Govaere locus_cov files lack {len(missing)} features; recompute them")
        ncov = pd.concat([ncov, g_ncov.reindex(ncov.index)], axis=1)
        genes = pd.concat([genes, g_genes], axis=1)
        qc = pd.concat([qc, g_qc])
        pheno = pd.concat([pheno, g_pheno])
    res, sf = chronic.d1_chronic(ncov, pheno,
                                 genes, qc, symbol_to_id, gtf.loc[gtf["gene_type"] == "protein_coding", "gene_id"],
                                 parents=parents, pseudocount=args.pseudocount)
    sf.rename("size_factor").to_csv(os.path.join(d, "size_factors.tsv"), sep="\t", index_label="sample_id")
    res = res.reindex(feats.index)
    path = os.path.join(args.out_dir, "D1.chronic.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; tested {res['D1_chronic'].notna().sum()}, FDR < 0.05: {(res['chronic_fdr'] < 0.05).sum()}, "
          f"parent-driven: {res['parent_driven'].sum()}")


def cmd_nmd(args):
    """M2 NMD sensitivity from the HeLa NMD knock-downs (quant_nmd/): manifest (needs network) or statistics."""
    d = os.path.join(args.out_dir, "quant_nmd")
    if args.manifest:
        os.makedirs(d, exist_ok=True)
        nmd.manifest().to_csv(os.path.join(d, "samples.tsv"), sep="\t", index=False)
        return
    ncov, samples, feats = expression.load(args.out_dir, "quant_nmd")
    ncov = ncov.loc[ncov.index.isin(expression.quantifiable(feats))]
    genes = quant.read_gene_sums(resources.fetch(f"sra_{nmd.PROJECT}_genes", args.data_dir,
                                                 url=quant.gene_sums_url("sra", nmd.PROJECT)), "sra", samples)
    pc = expression.protein_coding_ids(resources.fetch("gencode_gtf", args.data_dir))
    res = nmd.sensitivity(ncov, samples, genes[genes.index.isin(pc)], pseudocount=args.pseudocount)
    res = res.reindex(feats.index)
    path = os.path.join(args.out_dir, "nmd.loci.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; tested {res['nmd_log2fc'].notna().sum()}, up at FDR < 0.05: "
          f"{((res['nmd_fdr'] < 0.05) & (res['nmd_log2fc'] > 0)).sum()}")


def cmd_hcc(args):
    """D1 HCC replication in independent tumour / adjacent-normal cohorts (quant_hcc/): manifest (needs network) or
    statistics."""
    d = os.path.join(args.out_dir, "quant_hcc")
    if args.manifest:
        os.makedirs(d, exist_ok=True)
        man = hcc.manifest()
        man.to_csv(os.path.join(d, "samples.tsv"), sep="\t", index=False)
        print(man.groupby(["cohort", "tissue"])["patient"].nunique().to_string())
        return
    ncov, samples, feats = expression.load(args.out_dir, "quant_hcc")
    ncov = ncov.loc[ncov.index.isin(expression.quantifiable(feats))]
    pc = expression.protein_coding_ids(resources.fetch("gencode_gtf", args.data_dir))
    genes = pd.concat([quant.read_gene_sums(resources.fetch(f"sra_{p}_genes", args.data_dir,
                                                            url=quant.gene_sums_url("sra", p)), "sra", samples)
                       for p in hcc.COHORTS.values()], axis=1)
    sf = hcc.size_factors(genes[genes.index.isin(pc)], samples)
    sf.rename("size_factor").to_csv(os.path.join(d, "size_factors.tsv"), sep="\t", index_label="sample_id")
    res = hcc.replication(expression.scaled(ncov, sf), samples, pseudocount=args.pseudocount).reindex(feats.index)
    path = os.path.join(args.out_dir, "D1.hcc_replication.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; tested {res['hcc_rep_r'].notna().sum()}, FDR < 0.05: {(res['hcc_rep_fdr'] < 0.05).sum()}")


def cmd_literature(args):
    names = feature_names(args.results_dir, args.candidates)
    aliases = literature.hgnc_aliases(resources.fetch("hgnc", args.data_dir))
    res = literature.d5(names, aliases, os.path.join(args.data_dir, "europepmc_counts.json"), threads=args.threads)
    path = os.path.join(args.out_dir, "D5.literature.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}")


def cmd_genetics(args):
    feats = pd.read_csv(os.path.join(args.out_dir, "quant", "features.tsv"), sep="\t", index_col=0)
    assoc = genetics.read_liver_associations(resources.fetch("gwas_catalog", args.data_dir))
    res = genetics.d4(feats[["chrom", "start0", "end1"]], assoc, flank=args.flank)
    path = os.path.join(args.out_dir, "D4.gwas.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; {len(assoc)} liver-trait associations; features with a hit: {(res['gwas_n'] > 0).sum()}")

    # GTEx v8 liver eQTLs (reported, not scored); the two Liver files are extracted from the GTEx tar beforehand:
    # curl -s <eqtl.GTEX_TAR> | tar -x -C <data_dir>/gtex_eqtl --wildcards '*Liver.v8.egenes*' '*Liver.v8.signif*'
    g = os.path.join(args.data_dir, "gtex_eqtl")
    egenes = pd.read_csv(os.path.join(g, eqtl.EGENES), sep="\t")
    pairs = pd.read_csv(os.path.join(g, eqtl.PAIRS), sep="\t", usecols=["variant_id", "gene_id", "pval_nominal", "slope"])
    loci = inputs.loci(inputs.load(args.results_dir, args.candidates))
    res = eqtl.per_feature(eqtl.feature_genes(loci, controls.control_genes(), egenes), egenes, pairs, assoc)
    res.to_csv(os.path.join(args.out_dir, "D4.eqtl.tsv"), sep="\t", index_label="feature_id")
    print(f"eqtl: tested {len(res)}, eGenes {res['eqtl_egene'].sum()}, GWAS-linked {(res['eqtl_gwas_mlog10p'] > 0).sum()}")


def cmd_somatic(args):
    _, expr, samples, feats, _ = tcga_gtex_expression(args.out_dir, args.pseudocount)
    segments = pd.read_csv(resources.fetch("tcga_lihc_cn", args.data_dir), sep="\t")
    cn = somatic.cn_at(feats[["chrom", "start0", "end1"]], segments)
    purity = somatic.read_purity(resources.fetch("tcga_purity", args.data_dir))
    q = expression.quantifiable(feats)
    res = somatic.d3(expr.loc[expr.index.isin(q)], samples, cn, purity=purity).reindex(feats.index)
    path = os.path.join(args.out_dir, "D3.somatic.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; CN tested {res['cn_rho'].notna().sum()}")


def fetch_encode(data_dir):
    """ENCODE HepG2 fractionation bigWigs, eCLIP IDR peaks and the RBP localisation table,
    histone ChIP-seq bigWigs; returns the eCLIP file table."""
    d = os.path.join(data_dir, "encode")
    for files in encode.BIGWIGS.values():
        for acc in (a for rep in files for a in rep):
            resources.fetch(acc, d, url=encode.ENCODE_DOWNLOAD.format(acc))
    ed = os.path.join(d, "eclip")
    files = eclip.idr_files(json.load(open(resources.fetch("eclip_search", ed, url=eclip.ENCODE_SEARCH))))
    for acc in files["accession"]:
        resources.fetch(acc, ed, url=eclip.ENCODE_BED.format(acc))
    resources.fetch("rbp_localisation", ed, url=eclip.RBP_ANNOTATION)
    cd = os.path.join(d, "chromatin")
    for _, acc in chromatin.BIGWIGS.values():
        resources.fetch(acc, cd, url=chromatin.ENCODE_DOWNLOAD.format(acc))
    return files


def cmd_fetch(args):
    """Every public download used by the steps (needs network: run on the login node). Files already in
    --data-dir are kept; files without a stable URL are only checked (resources.MANUAL)."""
    for key in resources.URLS:
        resources.fetch(key, args.data_dir)
    for s, p in tcga_gtex_projects(args.data_dir):
        resources.fetch(f"{s}_{p}_genes", args.data_dir, url=quant.gene_sums_url(s, p))
    for p in hcc.COHORTS.values():
        resources.fetch(f"sra_{p}_genes", args.data_dir, url=quant.gene_sums_url("sra", p))
    for p in chronic.COHORTS.values():
        resources.fetch(f"sra_{p}_genes", args.data_dir, url=quant.gene_sums_url("sra", p))
        resources.fetch(f"sra_{p}_qc", args.data_dir, url=quant._metadata_url("sra", p, "recount_qc"))
    resources.fetch(f"sra_{nmd.PROJECT}_genes", args.data_dir, url=quant.gene_sums_url("sra", nmd.PROJECT))
    fetch_encode(args.data_dir)
    missing = resources.missing_manual(args.data_dir)
    for path, source in missing:
        print(f"MISSING {path}\n    get it from: {source}")
    if missing:
        raise SystemExit(f"{len(missing)} manually downloaded files are missing from {args.data_dir}")
    print(f"all resources present in {args.data_dir}")


def cmd_encode(args):
    d = os.path.join(args.data_dir, "encode")
    files = fetch_encode(args.data_dir)
    loci = inputs.loci(inputs.load(args.results_dir, args.candidates))
    gencode = quant.read_gtf_exons(resources.fetch("gencode_gtf", args.data_dir))
    # stranded data: only same-strand GENCODE exons compete for reads (quant/regions.tsv is for unstranded data)
    regions, strands = eclip.stranded_regions(
        quant.read_gtf_exons(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.gtf")), gencode, loci,
        controls.control_genes())
    res = encode.features(regions, strands, d)
    path = os.path.join(args.out_dir, "encode.loci.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; nuclear ratio for {res['nuclear_log2ratio'].notna().sum()} features")

    ed = os.path.join(d, "eclip")
    groups = eclip.rbp_groups(resources.fetch("rbp_localisation", ed, url=eclip.RBP_ANNOTATION))
    res = eclip.features(regions, strands, eclip.load_peak_sets(files, ed), groups)
    res.to_csv(os.path.join(args.out_dir, "eclip.loci.tsv"), sep="\t", index_label="feature_id")
    print(f"eclip: bound by >= 1 RBP: {(res['eclip_n_rbps'] > 0).sum()} of {len(res)} features")

    cd = os.path.join(d, "chromatin")
    tx =pd.read_csv(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.tsv"), sep="\t")
    tss = chromatin.dominant_tss(tx, resources.fetch("gencode_gtf", args.data_dir), controls.control_genes())
    tss = tss[tss.index.isin(loci.index) | tss.index.str.startswith("ctrl:")]
    res = chromatin.features(tss, cd)
    res.to_csv(os.path.join(args.out_dir, "chromatin.loci.tsv"), sep="\t", index_label="feature_id")
    print(f"chromatin: active TSS {(res['chrom_active'] > 0).sum()}, enhancer-like {(res['enhancer_score'] > 0).sum()}")


COEXPRESSION_PROJECTS = (("tcga", "LIHC"), ("gtex", "LIVER"))   # samples of the M3 co-expression


def neighbour_coexpression(args, loci_tx, gencode):
    """M3 co-expression of each quantifiable locus with its nearest protein-coding gene (TCGA-LIHC, GTEx liver)."""
    _, expr, samples, feats, sf = tcga_gtex_expression(args.out_dir, args.pseudocount)
    samples = samples[[(s, p) in COEXPRESSION_PROJECTS for s, p in zip(samples["source"], samples["project"])]]
    genes = pd.concat([quant.read_gene_sums(resources.fetch(f"{s}_{p}_genes", args.data_dir,
                                                            url=quant.gene_sums_url(s, p)), s, samples)
                       for s, p in COEXPRESSION_PROJECTS], axis=1)
    genes = expression.scaled(genes, sf)
    pc = gencode[gencode["gene_type"] == "protein_coding"].assign(gene_id=lambda d: d["gene_id"].str.split(".").str[0])
    pc_tss = mechanism.gene_spans(pc, "gene_id")
    loci_tss = mechanism.gene_spans(loci_tx, "gene_id")
    loci_tss = loci_tss[loci_tss.index.isin(expression.quantifiable(feats))]
    neighbours = mechanism.nearest_coding(loci_tss, pc_tss)
    null_pairs = mechanism.nearest_coding(pc_tss, pc_tss, exclude_self=True)
    cols = list(samples["sample_id"])
    return mechanism.coexpression(expr.loc[loci_tss.index, cols], genes[cols], neighbours, null_pairs)


def cmd_mechanism(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    orf_qc = pd.read_csv(os.path.join(args.out_dir, "qc.orfs.tsv"), sep="\t").set_index("orf_id").reindex(orfs["orf_id"])
    orf_qc.index = orfs.index
    peptides = pd.read_csv(os.path.join(args.out_dir, "peptides.orfs.tsv"), sep="\t", index_col=0)
    tx_exons = quant.read_gtf_exons(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.gtf"))
    loci_tx = pd.read_csv(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.tsv"), sep="\t")
    gencode = quant.read_gtf_exons(resources.fetch("gencode_gtf", args.data_dir))
    coexpr = neighbour_coexpression(args, loci_tx, gencode)
    coexpr.to_csv(os.path.join(args.out_dir, "coexpression.loci.tsv"), sep="\t", index_label="feature_id")
    feats = mechanism.features(orfs, orf_qc, peptides, mechanism.stop_to_last_junction(orfs, tx_exons),
                               mechanism.neighbour_position(loci_tx, gencode),
                               mechanism.small_rna_hosts(loci_tx, mechanism.small_rna_genes(
                                   gencode, resources.read_wgrna(args.data_dir))),
                               coexpr["coexpr_excess"],
                               *(pd.read_csv(os.path.join(args.out_dir, f"{t}.loci.tsv"), sep="\t", index_col=0)
                                 for t in ("encode", "eclip", "chromatin")),
                               peptide_feats=pd.read_csv(os.path.join(args.out_dir, "peptide_features.orfs.tsv"),
                                                         sep="\t", index_col=0),
                               orf_support=pd.read_csv(os.path.join(args.out_dir, "support.orfs.tsv"), sep="\t",
                                                       index_col=0),
                               nmd=pd.read_csv(os.path.join(args.out_dir, "nmd.loci.tsv"), sep="\t", index_col=0))
    res = mechanism.class_table(feats, pd.Series(True, index=orfs.index), orfs["locus_id"])
    out = pd.concat([orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate"]],
                     feats.add_prefix("raw_"), res], axis=1)
    path = os.path.join(args.out_dir, "mechanism.orfs.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"wrote {path}")
    print("dominant class, candidates:", out.loc[out["is_candidate"], "dominant_class"].value_counts().to_dict())


def cmd_disease(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    loci = inputs.loci(orfs)
    ncov, expr, samples, feats, sf = tcga_gtex_expression(args.out_dir, args.pseudocount)
    print(samples.groupby(["source", "project"]).size().to_string())
    cdr = pd.read_excel(resources.fetch("tcga_cdr", args.data_dir), sheet_name=0)

    q = expr.index.isin(expression.quantifiable(feats))
    d1 = disease.d1_hcc(expr.loc[q], samples)
    d2 = disease.d2_survival(expr.loc[q], samples, cdr, ncov=ncov.loc[q], size_factors=sf)
    cell = loci[["log2TPM_HepG2", "log2TPM_Huh7"]].astype(float).max(axis=1)
    d6 = integrate.liver_context(expr, samples, cell)
    d4 = pd.read_csv(os.path.join(args.out_dir, "D4.gwas.tsv"), sep="\t", index_col=0)
    d5 = pd.read_csv(os.path.join(args.out_dir, "D5.literature.tsv"), sep="\t", index_col=0)
    d3 = pd.read_csv(os.path.join(args.out_dir, "D3.somatic.tsv"), sep="\t", index_col=0)
    eqtl_table = pd.read_csv(os.path.join(args.out_dir, "D4.eqtl.tsv"), sep="\t", index_col=0)
    eqtl_table["D4_eqtl_supported"] = eqtl_table["eqtl_gwas_mlog10p"] > 0
    replication = pd.read_csv(os.path.join(args.out_dir, "D1.hcc_replication.tsv"), sep="\t", index_col=0)
    replication["hcc_replicated"] = hcc.replicated(replication, d1["D1_HCC"])

    stats = pd.DataFrame(index=feats.index)
    stats["D1_HCC"] = d1["D1_HCC"]
    d1c = pd.read_csv(os.path.join(args.out_dir, "D1.chronic.tsv"), sep="\t", index_col=0)
    stats["D1_chronic"] = d1c["D1_chronic"]
    stats["D2"] = d2["D2"]
    stats["D3"] = d3["D3"]
    stats["D4"] = d4["D4"]
    stats["D5"] = d5["D5"]
    for c in ("D4_hcc", "D4_chronic"):
        stats[c] = d4[c]
    for c in ("D5_hcc", "D5_chronic"):
        stats[c] = d5[c]
    stats = stats.join(d6)
    # loci that cannot be measured in unstranded public RNA-seq have no expression-based components
    not_quantifiable = ~feats["quantifiable"].reindex(stats.index).fillna(False).astype(bool)
    stats.loc[not_quantifiable, ["D1_HCC", "D1_chronic", "D2", "gtex_liver_median", "liver_specificity"]] = float("nan")
    is_bg = pd.Series(stats.index.isin(loci.index), index=stats.index)
    score = integrate.disease_score(stats, is_bg)

    tum_spec = integrate.tumour_specificity(expr, samples)
    tum_spec[not_quantifiable.reindex(tum_spec.index).fillna(False).astype(bool)] = float("nan")
    out = pd.concat([stats.add_prefix("stat_"), score, d1.drop(columns="D1_HCC"), d1c.drop(columns="D1_chronic"),
                     replication, d2.drop(columns="D2"),
                     d3.drop(columns="D3"), d4.drop(columns=["D4", "D4_hcc", "D4_chronic"]), eqtl_table,
                     d5.drop(columns=["D5", "D5_hcc", "D5_chronic", "symbol"]),
                     tum_spec.rename("tumour_specificity")], axis=1)
    out.insert(0, "symbol", loci["locus_name"].reindex(out.index).fillna(out.index.to_series().str.replace("ctrl:", "")))
    out.insert(1, "is_control", out.index.str.startswith("ctrl:"))
    out.insert(2, "is_candidate", loci["is_candidate"].reindex(out.index).fillna(False).astype(bool))
    out = out.sort_values("disease_score", ascending=False)
    path = os.path.join(args.out_dir, "disease.loci.tsv")
    out.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}")


def cmd_weights(args):
    path = os.path.join(args.out_dir, "disease.loci.tsv")
    dis = pd.read_csv(path, sep="\t", index_col=0)
    comp = dis[[f"comp_{c}" for c in weights.COMPONENTS]].rename(columns=lambda c: c[5:])
    y = weights.labels(dis[["symbol", "n_liver", "n_retracted"]], ~dis["is_control"], args.min_papers)
    w, summary = weights.fit(comp, y)
    weights.fitted_score(comp, w).rename("disease_score_fitted").to_csv(
        os.path.join(args.out_dir, "disease_score_fitted.tsv"), sep="\t", index_label="feature_id")
    out = pd.concat([w.rename("weight"), pd.Series(summary, name="value")])
    out.to_csv(os.path.join(args.out_dir, "weights.tsv"), sep="\t", header=["value"], index_label="key")
    print(out.to_string())


def cmd_report(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    read = lambda name, **kw: pd.read_csv(os.path.join(args.out_dir, name), sep="\t", **kw)
    orf_qc = read("qc.orfs.tsv").set_index("orf_id").reindex(orfs["orf_id"]).set_axis(orfs.index)
    mech = read("mechanism.orfs.tsv").set_index("orf_id").reindex(orfs["orf_id"]).set_axis(orfs.index)
    dis = read("disease.loci.tsv", index_col=0)
    if os.path.exists(os.path.join(args.out_dir, "disease_score_fitted.tsv")):
        dis = dis.drop(columns="disease_score_fitted", errors="ignore").join(
            read("disease_score_fitted.tsv", index_col=0))
    peptides = read("peptides.orfs.tsv", index_col=0)
    act = report.actionability(orfs, dis, peptides, read("hla.orfs.tsv", index_col=0))
    table = report.orf_table(orfs, orf_qc, dis, mech, act)
    table.to_csv(os.path.join(args.out_dir, "lincdisease.orfs.tsv"), sep="\t", index=False)
    cand = table[table["is_candidate"]].sort_values(["disease_score", "priority_score"], ascending=False)
    cand.to_csv(os.path.join(args.out_dir, "lincdisease.candidates.tsv"), sep="\t", index=False)
    screen = report.screen_supported(table, read("support.orfs.tsv"))
    screen.to_csv(os.path.join(args.out_dir, "lincdisease.screen_supported.tsv"), sep="\t", index=False)
    by_name = mech.assign(rank=orfs["rank"]).sort_values("rank").drop_duplicates("locus_name").set_index("locus_name")
    ctrl = report.control_report(dis, by_name, read("hla.control.tsv"))
    ctrl.to_csv(os.path.join(args.out_dir, "controls.tsv"), sep="\t", index=False)
    print(ctrl.to_string(index=False))
    print(f"wrote lincdisease.orfs.tsv, lincdisease.candidates.tsv, lincdisease.screen_supported.tsv, controls.tsv "
          f"to {args.out_dir}")


def cmd_html(args):
    """lincdisease.table.html: candidates and screen-supported ORFs, sortable, summary scores first."""
    read = lambda name, **kw: pd.read_csv(os.path.join(args.out_dir, name), sep="\t", low_memory=False, **kw)
    t = table_html.rows(read("lincdisease.orfs.tsv"), read("support.orfs.tsv"), read("mechanism.orfs.tsv"),
                        read("disease.loci.tsv", index_col=0), read("peptides.orfs.tsv", index_col=0),
                        read("hla.orfs.tsv", index_col=0))
    path = os.path.join(args.out_dir, "lincdisease.table.html")
    with open(path, "w") as fh:
        fh.write(table_html.build(t))
    print(f"wrote {path}: {len(t)} ORFs")


def main():
    p = argparse.ArgumentParser(prog="lincdisease", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("qc", help="Step 1: QC flags, copy-attribution confidence, confidence tier")
    _common(q)
    q.add_argument("--peptides-dir", required=True,
                   help="directory with proteomics/ and peptidomics/ results of the lincorfz notebook")
    q.set_defaults(func=cmd_qc)

    hm = sub.add_parser("homology", help="Step 1: DIAMOND homology to annotated proteins")
    _common(hm)
    hm.add_argument("--diamond", default="diamond")
    hm.add_argument("--threads", type=int, default=4)
    hm.set_defaults(func=cmd_homology)

    pp = sub.add_parser("peptide", help="Step 4 M1 sequence features and actionability HLA-I presentation")
    _common(pp)
    pp.add_argument("--tool-env", required=True, help="conda env with hmmer, metapredict and mhcflurry")
    pp.add_argument("--threads", type=int, default=4)
    pp.add_argument("--skip-hla", action="store_true")
    pp.set_defaults(func=cmd_peptide)

    su = sub.add_parser("support", help="Step 4 M1: ORF vs exon PhyloP contrast, published CRISPR screen hits")
    _common(su)
    su.add_argument("--phylop", default=peptide.PHYLOP_URL, help="PhyloP100way bigWig (path or URL)")
    su.add_argument("--threads", type=int, default=4)
    su.set_defaults(func=cmd_support)

    fe = sub.add_parser("fetch", help="All public downloads into --data-dir (needs network); checks hand-placed files")
    fe.add_argument("--data-dir", required=True)
    fe.set_defaults(func=cmd_fetch)

    r = sub.add_parser("regions", help="Step 2: quantifiable exonic regions of loci and controls")
    _common(r)
    r.add_argument("--min-quant-bases", type=int, default=100)
    r.set_defaults(func=cmd_regions)

    m = sub.add_parser("manifest", help="Step 2: recount3 samples (TCGA-LIHC, GTEx)")
    m.add_argument("--data-dir", required=True)
    m.add_argument("--gtex-per-tissue", type=int, default=10)
    m.add_argument("--seed", type=int, default=1)
    m.set_defaults(func=cmd_manifest)

    s = sub.add_parser("quant", help="Step 2: stream recount3 bigWigs and sum coverage over the regions")
    s.add_argument("--project", nargs="*", help="restrict to these recount3 projects")
    s.add_argument("--worker", type=int, default=0)
    s.add_argument("--n-workers", type=int, default=1)
    s.add_argument("--scratch", help="directory for the temporary bigWig downloads")
    s.add_argument("--sample-set", default="quant",
                   help="subdirectory of --out-dir with samples.tsv (quant: TCGA/GTEx; quant_sra: chronic cohorts)")
    s.add_argument("--features", nargs="*", help="quantify only these features (added later, e.g. a new control)")
    s.add_argument("--sums-dir", default="sums", help="subdirectory of quant/ for the sums (sums_<name> for added "
                                                       "features; joined to quant/sums/ by sample)")
    s.set_defaults(func=cmd_quant)

    sc = sub.add_parser("scaling", help="Step 2: size factors of the TCGA-LIHC / GTEx samples (recount3 gene sums)")
    sc.add_argument("--data-dir", required=True)
    sc.set_defaults(func=cmd_scaling)

    ch = sub.add_parser("chronic", help="Step 3, D1 chronic liver disease (SRA cohorts in quant_sra/)")
    ch.add_argument("--manifest", action="store_true", help="write quant_sra/samples.tsv (needs network)")
    ch.add_argument("--pseudocount", type=float, default=0.05)
    ch.add_argument("--data-dir", help="download cache (recount3 gene sums and QC, GENCODE); not needed with --manifest")
    ch.add_argument("--govaere-dir", help="Govaere directory (locus_cov/, counts/, star_logs/); its runs enter as "
                                          "the cohort 'govaere'")
    ch.add_argument("--govaere-exclude", nargs="*", default=[], help="Govaere runs left out (config govaere.exclude)")
    ch.add_argument("--govaere-runs", help="ENA run table of the Govaere runs")
    ch.add_argument("--govaere-soft", help="GEO sample annotation of GSE135251 (acc.cgi text format)")
    ch.set_defaults(func=cmd_chronic)

    nm = sub.add_parser("nmd", help="Step 4, M2: NMD sensitivity from HeLa NMD knock-downs (quant_nmd/)")
    nm.add_argument("--manifest", action="store_true", help="write quant_nmd/samples.tsv (needs network)")
    nm.add_argument("--pseudocount", type=float, default=0.05)
    nm.add_argument("--data-dir", help="download cache (recount3 gene sums, GENCODE); not needed with --manifest")
    nm.set_defaults(func=cmd_nmd)

    hc = sub.add_parser("hcc", help="Step 3, D1 HCC replication in independent cohorts (quant_hcc/)")
    hc.add_argument("--manifest", action="store_true", help="write quant_hcc/samples.tsv (needs network)")
    hc.add_argument("--pseudocount", type=float, default=0.05)
    hc.add_argument("--data-dir", help="download cache (recount3 gene sums, GENCODE); not needed with --manifest")
    hc.set_defaults(func=cmd_hcc)

    lt = sub.add_parser("literature", help="Step 3, D5: Europe PMC counts")
    _common(lt)
    lt.add_argument("--threads", type=int, default=4)
    lt.set_defaults(func=cmd_literature)

    gw = sub.add_parser("genetics", help="Step 3, D4: GWAS Catalog liver traits near each feature, GTEx liver eQTLs")
    _common(gw)
    gw.add_argument("--flank", type=int, default=50_000)
    gw.set_defaults(func=cmd_genetics)

    so = sub.add_parser("somatic", help="Step 3, D3: somatic copy number in TCGA-LIHC")
    _common(so)
    so.add_argument("--pseudocount", type=float, default=0.05, help="added to ncov before log2")
    so.set_defaults(func=cmd_somatic)

    en = sub.add_parser("encode", help="Step 4: ENCODE HepG2 fractionation, eCLIP and chromatin")
    _common(en)
    en.set_defaults(func=cmd_encode)

    mc = sub.add_parser("mechanism", help="Step 4: mechanism class features and scores")
    _common(mc)
    mc.add_argument("--pseudocount", type=float, default=0.05, help="added to ncov before log2")
    mc.set_defaults(func=cmd_mechanism)

    ds = sub.add_parser("disease", help="Step 3: disease components D1, D2, D6, collects D3-D5, locus disease score")
    _common(ds)
    ds.add_argument("--pseudocount", type=float, default=0.05, help="added to ncov before log2")
    ds.set_defaults(func=cmd_disease)

    wt = sub.add_parser("weights", help="Fitted disease-score weights (alternative score, reported)")
    _common(wt)
    wt.add_argument("--min-papers", type=int, default=weights.MIN_PAPERS,
                    help="liver papers (non-retracted) for an established liver-disease lncRNA")
    wt.set_defaults(func=cmd_weights)

    rp = sub.add_parser("report", help="Outputs: per-ORF table, candidates, control panel")
    _common(rp)
    rp.set_defaults(func=cmd_report)

    ht = sub.add_parser("html", help="Outputs: sortable HTML table of the candidates and screen-supported ORFs")
    ht.add_argument("--data-dir", help="unused; accepted for a uniform command line")
    ht.set_defaults(func=cmd_html)

    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main()
