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
import json
import os

import numpy as np
import pandas as pd

from . import (chromatin, chronic, controls, disease, eclip, encode, eqtl, expression, genetics, hcc, hla, homology,
               inputs, integrate, literature, mechanism, nmd, peptide, qc, quant, report, resources, screens, somatic,
               table_html, weights)


# ------------------------------------------------------------------ reading and writing

def read(out_dir, name, **kw):
    """A table written by an earlier step (tab-separated, in out_dir)."""
    return pd.read_csv(os.path.join(out_dir, name), sep="\t", **kw)


def write(df, out_dir, name, **kw):
    """Write a tab-separated table into out_dir; returns its path."""
    path = os.path.join(out_dir, name)
    df.to_csv(path, sep="\t", **kw)
    return path


def read_per_orf(out_dir, name, orfs):
    """A per-ORF table of an earlier step, in the row order (and index) of `orfs`: the step's row for each
    orfs["orf_id"] (NaN where it has none), relabelled with the index of `orfs`."""
    return read(out_dir, name).set_index("orf_id").reindex(orfs["orf_id"]).set_axis(orfs.index)


def lincrna_exons(results_dir):
    return quant.read_gtf_exons(os.path.join(results_dir, "lincRNAs", "lincRNA_transcripts.gtf"))


def lincrna_transcripts(results_dir):
    return pd.read_csv(os.path.join(results_dir, "lincRNAs", "lincRNA_transcripts.tsv"), sep="\t")


def gencode_exons(data_dir):
    return quant.read_gtf_exons(resources.fetch("gencode_gtf", data_dir))


def fetch_gene_sums(data_dir, source, project):
    """Path of the recount3 gene sums of a project (downloaded once)."""
    return resources.fetch(f"{source}_{project}_genes", data_dir, url=quant.gene_sums_url(source, project))


def gene_sums(data_dir, projects, samples):
    """AUC-scaled recount3 gene sums (genes x samples) of the given samples of the (source, project) pairs. Raises
    if a sample of these projects has none (it would drop out of the size factors and every test without notice)."""
    genes = pd.concat([quant.read_gene_sums(fetch_gene_sums(data_dir, s, p), s, samples) for s, p in projects],
                      axis=1)
    wanted = samples[[sp in set(projects) for sp in zip(samples["source"], samples["project"])]]["sample_id"]
    missing = wanted[~wanted.isin(genes.columns)]
    if len(missing):
        raise ValueError(f"no recount3 gene sums for {len(missing)} samples ({list(missing[:3])})")
    return genes


def quantifiable_only(df, feats):
    """Rows of the features measurable in unstranded public RNA-seq; the others are not tested (and stay out of the
    BH family)."""
    return df.loc[df.index.isin(expression.quantifiable(feats))]


def write_manifest(out_dir, sample_set, man):
    """<sample_set>/samples.tsv: the samples a `quant --sample-set` run streams."""
    os.makedirs(os.path.join(out_dir, sample_set), exist_ok=True)
    write(man, out_dir, os.path.join(sample_set, "samples.tsv"), index=False)


def write_size_factors(sf, out_dir, sample_set):
    write(sf.rename("size_factor"), out_dir, os.path.join(sample_set, "size_factors.tsv"), index_label="sample_id")


def feature_names(results_dir, candidates):
    """feature_id -> symbol for every locus (its locus_name) and control ('ctrl:<gene>')."""
    loci = inputs.loci(inputs.load(results_dir, candidates))
    names = loci["locus_name"].copy()
    for g in controls.control_genes():
        names[f"ctrl:{g}"] = g
    return names


# ------------------------------------------------------------------ Step 1: ORF quality

def read_peptides(pdir):
    """Specific and copy-specific peptides per ORF from the lincorfz proteomics and immunopeptidomics results."""
    def t(*parts):
        return pd.read_csv(os.path.join(pdir, *parts), sep="\t")
    return mechanism.specific_peptides(
        t("proteomics", "analysis", "lincORF_peptides.tsv"), t("proteomics", "fasta", "lincORF_entries.tsv"),
        t("peptidomics", "analysis", "immunopeptides.tsv"), t("peptidomics", "fasta", "lincORF_entries.tsv"))


def cmd_qc(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    peptides = read_peptides(args.peptides_dir)
    write(peptides, args.out_dir, "peptides.orfs.tsv", index_label="orf_id")
    kw = dict(segdups=resources.read_segdups(args.data_dir), false_dups=resources.read_false_dups(args.data_dir),
              copy_peptides=peptides["copy_specific"])
    # quantifiability and homology are optional (the workflow always has them)
    if os.path.exists(os.path.join(args.out_dir, "quant", "features.tsv")):
        feats = read(args.out_dir, "quant/features.tsv", index_col=0)
        kw["quantifiable"] = orfs.set_index("orf_id")["locus_id"].map(feats["quantifiable"])
    hom = None
    if os.path.exists(os.path.join(args.out_dir, "homology.orfs.tsv")):
        hom = read(args.out_dir, "homology.orfs.tsv")
        kw["homology"] = hom.loc[hom["parent"].notna(), ["orf_id", "parent", "pident", "qcov"]]
    res = qc.run_qc(orfs, **kw)
    if hom is not None:
        res["homolog_parent_coverage"] = hom.set_index("orf_id").reindex(orfs["orf_id"])["parent_coverage"].to_numpy()
    out = pd.concat([orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate"]],
                     res.drop(columns="orf_id")], axis=1)
    path = write(out, args.out_dir, "qc.orfs.tsv", index=False)
    print(f"wrote {path}: {len(out)} ORFs")
    print("tiers, all ORFs:", out["confidence_tier"].value_counts().to_dict())
    print("tiers, candidates:", out.loc[out["is_candidate"], "confidence_tier"].value_counts().to_dict())


def cmd_homology(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    res = homology.run(orfs, gencode_exons(args.data_dir), args.data_dir, os.path.join(args.out_dir, "homology"),
                       diamond=args.diamond, threads=args.threads)
    path = write(res, args.out_dir, "homology.orfs.tsv")
    print(f"wrote {path}: with homolog {res['parent'].notna().sum()}")


# ------------------------------------------------------------------ Step 2: quantification

TCGA_COLUMNS = ["tcga_barcode", "cgc_sample_sample_type", "gdc_cases.case_id", "gdc_cases.submitter_id"]
GTEX_COLUMNS = ["SUBJID", "SMTS", "SMTSD", "SEX", "AGE"]


def tcga_gtex_projects(data_dir):
    """(source, project) of the quant/ sample set: TCGA-LIHC, then every GTEx tissue."""
    projects = pd.read_csv(resources.fetch("gtex_projects", data_dir), sep="\t")["project"].unique()
    return [("tcga", "LIHC")] + [("gtex", p) for p in sorted(set(projects) - {"STUDY_NA"})]   # STUDY_NA: unknown


def cmd_regions(args):
    loci = inputs.loci(inputs.load(args.results_dir, args.candidates))
    regions, feats = quant.build_regions(lincrna_exons(args.results_dir), gencode_exons(args.data_dir), loci,
                                         controls.control_genes(), min_quant_bases=args.min_quant_bases)
    os.makedirs(os.path.join(args.out_dir, "quant"), exist_ok=True)
    write(regions, args.out_dir, "quant/regions.tsv", index=False)
    write(feats, args.out_dir, "quant/features.tsv", index_label="feature_id")
    print(f"{len(feats)} features, {len(regions)} regions; not quantifiable: {(~feats['quantifiable']).sum()}")


def cmd_manifest(args):
    """TCGA-LIHC (all samples), GTEx liver (all) and up to --gtex-per-tissue samples of every other GTEx tissue."""
    parts = [quant.manifest("tcga", "LIHC", TCGA_COLUMNS)]
    for _, project in tcga_gtex_projects(args.data_dir)[1:]:
        m = quant.manifest("gtex", project, GTEX_COLUMNS)
        if project != "LIVER" and len(m) > args.gtex_per_tissue:
            m = m.sample(args.gtex_per_tissue, random_state=args.seed)
        parts.append(m)
    man = pd.concat(parts, ignore_index=True)
    write_manifest(args.out_dir, "quant", man)
    print(man.groupby(["source", "project"]).size().to_string())


def cmd_quant(args):
    regions = read(args.out_dir, "quant/regions.tsv")
    d = os.path.join(args.out_dir, args.sample_set)
    samples = read(d, "samples.tsv", low_memory=False)
    if args.project:
        samples = samples[samples["project"].isin(args.project)]
    if args.features:
        regions = regions[regions["feature_id"].isin(args.features)]
    samples = samples.iloc[args.worker::args.n_workers]
    path = quant.run(samples, regions, os.path.join(d, args.sums_dir), worker=f"{args.worker}", scratch=args.scratch)
    print(f"done: {path}")


def cmd_scaling(args):
    """Size factors of the TCGA-LIHC and GTEx samples (quant/size_factors.tsv), over all of them together, since
    tumours, adjacent normals and the GTEx tissues are compared with each other (D1, D6, tumour specificity)."""
    _, samples, _ = expression.load(args.out_dir)          # the samples with coverage sums
    pc = expression.protein_coding_ids(resources.fetch("gencode_gtf", args.data_dir))
    present = set(zip(samples["source"], samples["project"]))
    genes = gene_sums(args.data_dir, [sp for sp in tcga_gtex_projects(args.data_dir) if sp in present], samples)
    sf = expression.size_factors(genes[genes.index.isin(pc)])
    write_size_factors(sf, args.out_dir, "quant")
    print(f"wrote {os.path.join(args.out_dir, 'quant', 'size_factors.tsv')}; median by group:")
    s = samples.set_index("sample_id")
    group = np.where(s["source"] == "tcga", "TCGA " + s["cgc_sample_sample_type"].astype(str), "GTEx " + s["project"])
    print(sf.groupby(pd.Series(group, index=s.index).reindex(sf.index)).median().round(2).to_string())


def tcga_gtex_expression(out_dir, pseudocount):
    """Size-factor scaled ncov of the TCGA-LIHC / GTEx samples, log2(ncov + pseudocount), samples, features and
    the size factors."""
    ncov, samples, feats = expression.load(out_dir)
    sf = expression.read_size_factors(out_dir)
    ncov = expression.scaled(ncov, sf)
    samples = samples[samples["sample_id"].isin(ncov.columns)]
    return ncov, disease.log_expr(ncov, pseudocount), samples, feats, sf


# ------------------------------------------------------------------ Step 3: disease components

def parent_genes(out_dir, symbol_to_id):
    """locus_id -> gene id of its parent gene: the homology parent of its ORF with the highest identity."""
    orf_qc = read(out_dir, "qc.orfs.tsv").dropna(subset=["homolog_parent"])
    best = orf_qc.sort_values("homolog_pident", ascending=False).drop_duplicates("locus_id").set_index("locus_id")
    return best["homolog_parent"].map(symbol_to_id).dropna()


def cmd_chronic(args):
    """D1 chronic liver disease from the SRA cohorts and the Govaere runs: manifest (needs network) or statistics
    on the quantified sums."""
    if args.manifest:
        write_manifest(args.out_dir, "quant_sra", chronic.manifest())
        return
    ncov, samples, feats = expression.load(args.out_dir, "quant_sra")
    ncov = quantifiable_only(ncov, feats)
    projects = chronic.COHORTS.values()
    genes = gene_sums(args.data_dir, [("sra", p) for p in projects], samples)
    run_qc = pd.concat([pd.read_csv(resources.fetch(f"sra_{p}_qc", args.data_dir,
                                                    url=quant.metadata_url("sra", p, "recount_qc")),
                                    sep="\t", low_memory=False) for p in projects])
    run_qc.index = "sra:" + run_qc["external_id"].astype(str)
    gtf = gencode_exons(args.data_dir).drop_duplicates("gene_id")
    gtf["gene_id"] = gtf["gene_id"].str.split(".").str[0]
    symbol_to_id = dict(zip(gtf["gene_name"], gtf["gene_id"]))
    pheno = samples.set_index("sample_id")[["cohort", "fibrosis_stage", "nas", "group"]]
    if args.govaere_dir:
        regions = read(args.out_dir, "quant/regions.tsv")
        g_ncov, g_genes, g_qc, g_pheno = chronic.govaere(args.govaere_dir, args.govaere_runs, args.govaere_soft,
                                                         regions, exclude=args.govaere_exclude)
        print(f"Govaere: {len(g_pheno)} runs")
        missing = ncov.index.difference(g_ncov.index)
        if len(missing):
            raise ValueError(f"Govaere locus_cov files lack {len(missing)} features; recompute them")
        ncov = pd.concat([ncov, g_ncov.reindex(ncov.index)], axis=1)
        genes = pd.concat([genes, g_genes], axis=1)
        run_qc = pd.concat([run_qc, g_qc])
        pheno = pd.concat([pheno, g_pheno])
    res, sf = chronic.d1_chronic(ncov, pheno, genes, run_qc, symbol_to_id,
                                 gtf.loc[gtf["gene_type"] == "protein_coding", "gene_id"],
                                 parents=parent_genes(args.out_dir, symbol_to_id), pseudocount=args.pseudocount)
    write_size_factors(sf, args.out_dir, "quant_sra")
    res = res.reindex(feats.index)
    path = write(res, args.out_dir, "D1.chronic.tsv", index_label="feature_id")
    print(f"wrote {path}; tested {res['D1_chronic'].notna().sum()}, FDR < 0.05: {(res['chronic_fdr'] < 0.05).sum()}, "
          f"parent-driven: {res['parent_driven'].sum()}")


def cmd_hcc(args):
    """D1 HCC replication in independent tumour / adjacent-normal cohorts (quant_hcc/): manifest (needs network) or
    statistics."""
    if args.manifest:
        man = hcc.manifest()
        write_manifest(args.out_dir, "quant_hcc", man)
        print(man.groupby(["cohort", "tissue"])["patient"].nunique().to_string())
        return
    ncov, samples, feats = expression.load(args.out_dir, "quant_hcc")
    ncov = quantifiable_only(ncov, feats)
    pc = expression.protein_coding_ids(resources.fetch("gencode_gtf", args.data_dir))
    genes = gene_sums(args.data_dir, [("sra", p) for p in hcc.COHORTS.values()], samples)
    sf = hcc.size_factors(genes[genes.index.isin(pc)], samples)
    write_size_factors(sf, args.out_dir, "quant_hcc")
    res = hcc.replication(expression.scaled(ncov, sf), samples, pseudocount=args.pseudocount).reindex(feats.index)
    path = write(res, args.out_dir, "D1.hcc_replication.tsv", index_label="feature_id")
    print(f"wrote {path}; tested {res['hcc_rep_r'].notna().sum()}, FDR < 0.05: {(res['hcc_rep_fdr'] < 0.05).sum()}")


def cmd_somatic(args):
    """D3: copy number vs expression of the loci, relative to the protein-coding genes around them."""
    _, expr, samples, feats, sf = tcga_gtex_expression(args.out_dir, args.pseudocount)
    segments = pd.read_csv(resources.fetch("tcga_lihc_cn", args.data_dir), sep="\t")
    cn = somatic.cn_at(feats[["chrom", "start0", "end1"]], segments)
    purity = somatic.read_purity(resources.fetch("tcga_purity", args.data_dir))
    # the neighbourhood: protein-coding genes of the same tumours
    tcga = samples[samples["source"] == "tcga"]
    pc = expression.protein_coding_ids(resources.fetch("gencode_gtf", args.data_dir))
    genes = gene_sums(args.data_dir, [("tcga", "LIHC")], tcga)
    genes = expression.scaled(genes[genes.index.isin(pc)], sf)
    spans = somatic.gene_spans(gencode_exons(args.data_dir), genes.index)
    genes = genes.loc[genes.index.intersection(spans.index)]
    rho = somatic.dosage_rho(np.log2(genes + 1), genes, tcga, somatic.cn_at(spans.loc[genes.index], segments),
                             purity=purity)
    neighbours = somatic.neighbour_median(feats[["chrom", "start0", "end1"]], spans, rho)
    print(f"protein-coding genes with a dosage rho: {rho.notna().sum()} (median {rho.median():.3f})")
    res = somatic.d3(quantifiable_only(expr, feats), samples, cn, purity=purity,
                     neighbour_rho=neighbours).reindex(feats.index)
    path = write(res, args.out_dir, "D3.somatic.tsv", index_label="feature_id")
    print(f"wrote {path}; CN tested {res['cn_rho'].notna().sum()}")


def cmd_genetics(args):
    feats = read(args.out_dir, "quant/features.tsv", index_col=0)
    assoc = genetics.read_liver_associations(resources.fetch("gwas_catalog", args.data_dir))
    res = genetics.d4(feats[["chrom", "start0", "end1"]], assoc, flank=args.flank)
    path = write(res, args.out_dir, "D4.gwas.tsv", index_label="feature_id")
    print(f"wrote {path}; {len(assoc)} liver-trait associations; features with a hit: {(res['gwas_n'] > 0).sum()}")

    # GTEx v8 liver eQTLs (reported, not scored); the two Liver files are hand-placed (resources.MANUAL)
    g = os.path.join(args.data_dir, "gtex_eqtl")
    egenes = pd.read_csv(os.path.join(g, eqtl.EGENES), sep="\t")
    pairs = pd.read_csv(os.path.join(g, eqtl.PAIRS), sep="\t",
                        usecols=["variant_id", "gene_id", "pval_nominal", "slope"])
    loci = inputs.loci(inputs.load(args.results_dir, args.candidates))
    res = eqtl.per_feature(eqtl.feature_genes(loci, controls.control_genes(), egenes), egenes, pairs, assoc)
    write(res, args.out_dir, "D4.eqtl.tsv", index_label="feature_id")
    print(f"eqtl: tested {len(res)}, eGenes {res['eqtl_egene'].sum()}, "
          f"GWAS-linked {(res['eqtl_gwas_mlog10p'] > 0).sum()}")


def cmd_literature(args):
    names = feature_names(args.results_dir, args.candidates)
    aliases = literature.hgnc_aliases(resources.fetch("hgnc", args.data_dir))
    res = literature.d5(names, aliases, os.path.join(args.data_dir, "europepmc_counts.json"), threads=args.threads)
    path = write(res, args.out_dir, "D5.literature.tsv", index_label="feature_id")
    print(f"wrote {path}")


def cmd_disease(args):
    """D1 HCC, D2 and D6 from the TCGA-LIHC / GTEx expression, the components of the other steps, the locus disease
    score (disease.loci.tsv)."""
    loci = inputs.loci(inputs.load(args.results_dir, args.candidates))
    ncov, expr, samples, feats, sf = tcga_gtex_expression(args.out_dir, args.pseudocount)
    print(samples.groupby(["source", "project"]).size().to_string())
    cdr = pd.read_excel(resources.fetch("tcga_cdr", args.data_dir), sheet_name=0)

    q = expr.index.isin(expression.quantifiable(feats))
    d1 = disease.d1_hcc(expr.loc[q], samples, ncov=ncov.loc[q])
    d2 = disease.d2_survival(expr.loc[q], samples, cdr, ncov=ncov.loc[q], size_factors=sf)
    cell = loci[["log2TPM_HepG2", "log2TPM_Huh7"]].astype(float).max(axis=1)
    d6 = integrate.liver_context(expr, samples, cell)
    d1c = read(args.out_dir, "D1.chronic.tsv", index_col=0)
    d3 = read(args.out_dir, "D3.somatic.tsv", index_col=0)
    d4 = read(args.out_dir, "D4.gwas.tsv", index_col=0)
    d5 = read(args.out_dir, "D5.literature.tsv", index_col=0)
    eqtl_table = read(args.out_dir, "D4.eqtl.tsv", index_col=0)
    eqtl_table["D4_eqtl_supported"] = eqtl_table["eqtl_gwas_mlog10p"] > 0
    replication = read(args.out_dir, "D1.hcc_replication.tsv", index_col=0)

    # the component statistics (stat_* in the output); D4_/D5_hcc and _chronic are the disease-specific parts
    stats = pd.concat([d1["D1_HCC"], d1c["D1_chronic"], d2["D2"], d3["D3"], d4["D4"], d5["D5"],
                       d4[["D4_hcc", "D4_chronic"]], d5[["D5_hcc", "D5_chronic"]]], axis=1).reindex(feats.index)
    stats = stats.join(d6)
    # loci that cannot be measured in unstranded public RNA-seq have no expression-based components
    not_quantifiable = ~feats["quantifiable"].reindex(stats.index).fillna(False).astype(bool)
    stats.loc[not_quantifiable, ["D1_HCC", "D1_chronic", "D2", "gtex_liver_median", "liver_specificity"]] = float("nan")
    is_bg = pd.Series(stats.index.isin(loci.index), index=stats.index)
    score = integrate.disease_score(stats, is_bg)
    replication["hcc_replicated"] = hcc.replicated(replication, d1["D1_HCC"], is_bg)

    tum_spec = integrate.tumour_specificity(expr, samples)
    tum_spec[not_quantifiable.reindex(tum_spec.index).fillna(False).astype(bool)] = float("nan")
    out = pd.concat([stats.add_prefix("stat_"), score, d1.drop(columns="D1_HCC"), d1c.drop(columns="D1_chronic"),
                     replication, d2.drop(columns="D2"),
                     d3.drop(columns="D3"), d4.drop(columns=["D4", "D4_hcc", "D4_chronic"]), eqtl_table,
                     d5.drop(columns=["D5", "D5_hcc", "D5_chronic", "symbol"]),
                     tum_spec.rename("tumour_specificity")], axis=1)
    gene = out.index.to_series().str.replace("ctrl:", "")
    out.insert(0, "symbol", loci["locus_name"].reindex(out.index).fillna(gene))
    out.insert(1, "is_control", out.index.str.startswith("ctrl:"))
    out.insert(2, "is_candidate", loci["is_candidate"].reindex(out.index).fillna(False).astype(bool))
    out = out.sort_values("disease_score", ascending=False)
    path = write(out, args.out_dir, "disease.loci.tsv", index_label="feature_id")
    print(f"wrote {path}")


def cmd_weights(args):
    dis = read(args.out_dir, "disease.loci.tsv", index_col=0)
    comp = dis[[f"comp_{c}" for c in weights.COMPONENTS]].rename(columns=lambda c: c[len("comp_"):])
    y = weights.labels(dis[["symbol", "n_liver", "n_retracted"]], ~dis["is_control"], args.min_papers)
    w, summary = weights.fit(comp, y)
    write(weights.fitted_score(comp, w).rename("disease_score_fitted"), args.out_dir, "disease_score_fitted.tsv",
          index_label="feature_id")
    out = pd.concat([w.rename("weight"), pd.Series(summary, name="value")])
    write(out, args.out_dir, "weights.tsv", header=["value"], index_label="key")
    print(out.to_string())


# ------------------------------------------------------------------ Step 4: mechanism

def cmd_peptide(args):
    """M1 sequence features (TM, signal-like, Pfam, disorder), GENCODE Ribo-seq ORF matches, MHCflurry HLA-I
    presentation. hmmsearch, metapredict and MHCflurry run in a separate environment (--tool-env)."""
    orfs = inputs.load(args.results_dir, args.candidates)
    env = os.path.join(args.tool_env, "bin")
    out = orfs[["orf_id"]].join(peptide.sequence_features(orfs))
    rd = os.path.join(args.data_dir, "riboseq_orfs")
    phase2 = peptide.read_bed12(os.path.join(rd, "phase2", "Ribo-seq_ORFs.comprehensive.bed"))
    phase1 = peptide.read_bed12(os.path.join(rd, "phase1", "Ribo-seq_ORFs.peptideins.bed"))
    out["gencode_riboseq_orfs"] = peptide.same_stop(orfs, phase2)
    out["gencode_riboseq_orf"] = out["gencode_riboseq_orfs"] != ""
    out["gencode_riboseq_peptide_support"] = peptide.same_stop(orfs, phase1) != ""
    seqs = orfs.set_index("orf_id")["aa_seq"].fillna("").str.rstrip("*")
    pf = peptide.pfam(seqs, os.path.join(env, "hmmsearch"), os.path.join(args.data_dir, "Pfam-A.hmm"), cpu=args.threads)
    out = out.join(pf, on="orf_id")
    out["pfam_n_domains"] = out["pfam_n_domains"].fillna(0).astype(int)
    out["pfam_domains"] = out["pfam_domains"].fillna("")
    out = out.join(peptide.disorder(seqs, os.path.join(env, "python")).add_prefix("pep_"), on="orf_id")
    write(out, args.out_dir, "peptide_features.orfs.tsv", index=False)
    print(f"peptide features: Ribo-seq catalog {out['gencode_riboseq_orf'].sum()}, "
          f"Pfam {(out['pfam_n_domains'] > 0).sum()}")
    if args.skip_hla:
        return
    # self peptides: those also found in the annotated proteomes (as for the homology search)
    hd = os.path.join(args.data_dir, "homology")
    proteins = [os.path.join(hd, f) for f in ("gencode.v42.pc_translations.fa.gz", "UP000005640_9606.fasta.gz",
                                              "UP000005640_9606_additional.fasta.gz")]
    res = hla.run(orfs, os.path.join(env, "python"), os.path.join(args.out_dir, "hla.orfs.tsv"),
                  env={"MHCFLURRY_DATA_DIR": os.path.join(args.data_dir, "mhcflurry")}, protein_fastas=proteins)
    print(f"HLA: ORFs with >= 1 non-self binder {(res['hla_any_n_binders'] > 0).sum()} of {len(res)}; "
          f"self-peptide binders left out: {res['hla_n_self_binders'].sum()}")


def cmd_support(args):
    """M1 conservation contrast (PhyloP100way, ORF vs the rest of its transcript's exons) and published CRISPR
    screen hits (Prensner et al. 2021 supplementary tables in <data_dir>/crispr_screens)."""
    orfs = inputs.load(args.results_dir, args.candidates)
    regions = read(args.out_dir, "quant/regions.tsv")
    cons = peptide.conservation_contrast(orfs, lincrna_exons(args.results_dir), bigwig=args.phylop,
                                         threads=args.threads, regions=regions)
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
    write(out, args.out_dir, "support.orfs.tsv", index=False)
    print(f"support: PhyloP contrast for {out['phylop_orf_minus_nonorf'].notna().sum()} ORFs; screen matches "
          f"Prensner {(out['prensner_orfs'] != '').sum()}, Chen {(out['chen_orfs'] != '').sum()}; combined values "
          f"{out['crispr_screen_value'].value_counts().to_dict()}")


def fetch_encode(data_dir):
    """ENCODE HepG2 fractionation bigWigs, eCLIP IDR peaks and the RBP localisation table,
    histone ChIP-seq bigWigs; returns the eCLIP file table."""
    d = os.path.join(data_dir, "encode")
    for files in encode.BIGWIGS.values():
        for acc in (a for rep in files for a in rep):
            resources.fetch(acc, d, url=encode.ENCODE_DOWNLOAD.format(acc))
    ed = os.path.join(d, "eclip")
    with open(resources.fetch("eclip_search", ed, url=eclip.ENCODE_SEARCH)) as fh:
        files = eclip.idr_files(json.load(fh))
    for acc in files["accession"]:
        resources.fetch(acc, ed, url=eclip.ENCODE_BED.format(acc))
    resources.fetch("rbp_localisation", ed, url=eclip.RBP_ANNOTATION)
    cd = os.path.join(d, "chromatin")
    for _, acc in chromatin.BIGWIGS.values():
        resources.fetch(acc, cd, url=chromatin.ENCODE_DOWNLOAD.format(acc))
    return files


def cmd_encode(args):
    """From ENCODE HepG2: nuclear / cytoplasmic localisation (fractionation; M5), eCLIP binding by nuclear (M4) and
    cytoplasmic (M5) RBPs, and chromatin at the TSS (M3 enhancer-like chromatin)."""
    d = os.path.join(args.data_dir, "encode")
    files = fetch_encode(args.data_dir)
    loci = inputs.loci(inputs.load(args.results_dir, args.candidates))
    # stranded data: only same-strand GENCODE exons compete for reads (quant/regions.tsv is for unstranded data)
    regions, strands = eclip.stranded_regions(lincrna_exons(args.results_dir), gencode_exons(args.data_dir), loci,
                                              controls.control_genes())
    res = encode.features(regions, strands, d)
    path = write(res, args.out_dir, "encode.loci.tsv", index_label="feature_id")
    print(f"wrote {path}; nuclear ratio for {res['nuclear_log2ratio'].notna().sum()} features")

    ed = os.path.join(d, "eclip")
    groups = eclip.rbp_groups(resources.fetch("rbp_localisation", ed, url=eclip.RBP_ANNOTATION))
    res = eclip.features(regions, strands, eclip.load_peak_sets(files, ed), groups)
    write(res, args.out_dir, "eclip.loci.tsv", index_label="feature_id")
    print(f"eclip: bound by >= 1 RBP: {(res['eclip_n_rbps'] > 0).sum()} of {len(res)} features")

    tss = chromatin.dominant_tss(lincrna_transcripts(args.results_dir), resources.fetch("gencode_gtf", args.data_dir),
                                 controls.control_genes())
    tss = tss[tss.index.isin(loci.index) | tss.index.str.startswith("ctrl:")]
    res = chromatin.features(tss, os.path.join(d, "chromatin"))
    write(res, args.out_dir, "chromatin.loci.tsv", index_label="feature_id")
    print(f"chromatin: active TSS {(res['chrom_active'] > 0).sum()}, enhancer-like {(res['enhancer_score'] > 0).sum()}")


def cmd_nmd(args):
    """M2 NMD sensitivity from the HeLa NMD knock-downs (quant_nmd/): manifest (needs network) or statistics."""
    if args.manifest:
        write_manifest(args.out_dir, "quant_nmd", nmd.manifest())
        return
    ncov, samples, feats = expression.load(args.out_dir, "quant_nmd")
    ncov = quantifiable_only(ncov, feats)
    genes = gene_sums(args.data_dir, [("sra", nmd.PROJECT)], samples)
    pc = expression.protein_coding_ids(resources.fetch("gencode_gtf", args.data_dir))
    res = nmd.sensitivity(ncov, samples, genes[genes.index.isin(pc)], pseudocount=args.pseudocount)
    res = res.reindex(feats.index)
    path = write(res, args.out_dir, "nmd.loci.tsv", index_label="feature_id")
    print(f"wrote {path}; tested {res['nmd_log2fc'].notna().sum()}, up at FDR < 0.05: "
          f"{((res['nmd_fdr'] < 0.05) & (res['nmd_log2fc'] > 0)).sum()}")


COEXPRESSION_PROJECTS = (("tcga", "LIHC"), ("gtex", "LIVER"))   # samples of the M3 co-expression


def neighbour_coexpression(out_dir, data_dir, pseudocount, loci_tx, gencode):
    """M3 co-expression of each quantifiable locus with its nearest protein-coding gene (TCGA-LIHC, GTEx liver)."""
    _, expr, samples, feats, sf = tcga_gtex_expression(out_dir, pseudocount)
    samples = samples[[sp in COEXPRESSION_PROJECTS for sp in zip(samples["source"], samples["project"])]]
    genes = expression.scaled(gene_sums(data_dir, COEXPRESSION_PROJECTS, samples), sf)
    pc = gencode[gencode["gene_type"] == "protein_coding"].assign(gene_id=lambda d: d["gene_id"].str.split(".").str[0])
    pc_tss = mechanism.gene_spans(pc, "gene_id")
    loci_tss = quantifiable_only(mechanism.gene_spans(loci_tx, "gene_id"), feats)
    neighbours = mechanism.nearest_coding(loci_tss, pc_tss)
    null_pairs = mechanism.nearest_coding(pc_tss, pc_tss, exclude_self=True)
    cols = list(samples["sample_id"])
    return mechanism.coexpression(expr.loc[loci_tss.index, cols], genes[cols], neighbours, null_pairs)


def cmd_mechanism(args):
    """Mechanism class features (M1-M6) of every ORF and the class scores (mechanism.orfs.tsv)."""
    orfs = inputs.load(args.results_dir, args.candidates)

    def indexed(name):
        return read(args.out_dir, name, index_col=0)

    loci_tx = lincrna_transcripts(args.results_dir)
    gencode = gencode_exons(args.data_dir)
    coexpr = neighbour_coexpression(args.out_dir, args.data_dir, args.pseudocount, loci_tx, gencode)
    write(coexpr, args.out_dir, "coexpression.loci.tsv", index_label="feature_id")
    small_rnas = mechanism.small_rna_genes(gencode, resources.read_wgrna(args.data_dir))
    feats = mechanism.features(orfs, read_per_orf(args.out_dir, "qc.orfs.tsv", orfs), indexed("peptides.orfs.tsv"),
                               mechanism.stop_to_last_junction(orfs, lincrna_exons(args.results_dir)),
                               mechanism.neighbour_position(loci_tx, gencode),
                               mechanism.small_rna_hosts(loci_tx, small_rnas),
                               coexpr["coexpr_excess"],
                               indexed("encode.loci.tsv"), indexed("eclip.loci.tsv"), indexed("chromatin.loci.tsv"),
                               peptide_feats=indexed("peptide_features.orfs.tsv"),
                               orf_support=indexed("support.orfs.tsv"), nmd=indexed("nmd.loci.tsv"))
    res = mechanism.class_table(feats, pd.Series(True, index=orfs.index), orfs["locus_id"])
    out = pd.concat([orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate"]],
                     feats.add_prefix("raw_"), res], axis=1)
    path = write(out, args.out_dir, "mechanism.orfs.tsv", index=False)
    print(f"wrote {path}")
    print("dominant class, candidates:", out.loc[out["is_candidate"], "dominant_class"].value_counts().to_dict())


# ------------------------------------------------------------------ outputs

def cmd_report(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    orf_qc = read_per_orf(args.out_dir, "qc.orfs.tsv", orfs)
    mech = read_per_orf(args.out_dir, "mechanism.orfs.tsv", orfs)
    dis = read(args.out_dir, "disease.loci.tsv", index_col=0)
    if os.path.exists(os.path.join(args.out_dir, "disease_score_fitted.tsv")):
        dis = dis.drop(columns="disease_score_fitted", errors="ignore").join(
            read(args.out_dir, "disease_score_fitted.tsv", index_col=0))
    peptides = read(args.out_dir, "peptides.orfs.tsv", index_col=0)
    act = report.actionability(orfs, dis, peptides, read(args.out_dir, "hla.orfs.tsv", index_col=0))
    table = report.orf_table(orfs, orf_qc, dis, mech, act)
    write(table, args.out_dir, "lincdisease.orfs.tsv", index=False)
    # headline: the disease score without prior literature (D5), unbiased by how well studied a locus is
    cand = table[table["is_candidate"]].sort_values(["disease_score_noD5", "priority_score"], ascending=False)
    write(cand, args.out_dir, "lincdisease.candidates.tsv", index=False)
    screen = report.screen_supported(table, read(args.out_dir, "support.orfs.tsv"))
    write(screen, args.out_dir, "lincdisease.screen_supported.tsv", index=False)
    # controls are matched by gene symbol, to the best-ranked ORF of the locus of that name
    by_name = mech.assign(rank=orfs["rank"]).sort_values("rank").drop_duplicates("locus_name").set_index("locus_name")
    ctrl = report.control_report(dis, by_name, read(args.out_dir, "hla.control.tsv"))
    write(ctrl, args.out_dir, "controls.tsv", index=False)
    print(ctrl.to_string(index=False))
    print(f"wrote lincdisease.orfs.tsv, lincdisease.candidates.tsv, lincdisease.screen_supported.tsv, controls.tsv "
          f"to {args.out_dir}")


def cmd_html(args):
    """lincdisease.table.html: candidates and screen-supported ORFs, sortable, summary scores first."""
    def t(name, **kw):
        return read(args.out_dir, name, low_memory=False, **kw)
    rows = table_html.rows(t("lincdisease.orfs.tsv"), t("support.orfs.tsv"), t("mechanism.orfs.tsv"),
                           t("disease.loci.tsv", index_col=0), t("peptides.orfs.tsv", index_col=0),
                           t("hla.orfs.tsv", index_col=0))
    path = os.path.join(args.out_dir, "lincdisease.table.html")
    with open(path, "w") as fh:
        fh.write(table_html.build(rows))
    print(f"wrote {path}: {len(rows)} ORFs")


def cmd_fetch(args):
    """Every public download used by the steps (needs network: run on the login node). Files already in
    --data-dir are kept; files without a stable URL are only checked (resources.MANUAL)."""
    for key in resources.URLS:
        resources.fetch(key, args.data_dir)
    sra = [*hcc.COHORTS.values(), *chronic.COHORTS.values(), nmd.PROJECT]
    for s, p in tcga_gtex_projects(args.data_dir) + [("sra", p) for p in sra]:
        fetch_gene_sums(args.data_dir, s, p)
    for p in chronic.COHORTS.values():
        resources.fetch(f"sra_{p}_qc", args.data_dir, url=quant.metadata_url("sra", p, "recount_qc"))
    fetch_encode(args.data_dir)
    missing = resources.missing_manual(args.data_dir)
    for path, source in missing:
        print(f"MISSING {path}\n    get it from: {source}")
    if missing:
        raise SystemExit(f"{len(missing)} manually downloaded files are missing from {args.data_dir}")
    print(f"all resources present in {args.data_dir}")


# ------------------------------------------------------------------ command line

def _run_inputs(p):
    """Options of the steps that start from the lincorfz run."""
    p.add_argument("--results-dir", required=True, help="lincorfz results/ directory (background)")
    p.add_argument("--candidates", required=True, help="candidate ORF table (lincorfz format)")
    p.add_argument("--data-dir", required=True, help="download cache for public resources")


def _threads(p):
    p.add_argument("--threads", type=int, default=4)


def _pseudocount(p):
    p.add_argument("--pseudocount", type=float, default=0.05, help="added to ncov before log2")


def _cohort_step(sub, name, sample_set, help):
    """chronic, nmd, hcc: --manifest writes <sample_set>/samples.tsv, otherwise the statistics."""
    p = sub.add_parser(name, help=help)
    p.add_argument("--manifest", action="store_true", help=f"write {sample_set}/samples.tsv (needs network)")
    _pseudocount(p)
    p.add_argument("--data-dir", help="download cache (recount3 gene sums, GENCODE); not needed with --manifest")
    return p


def main():
    p = argparse.ArgumentParser(prog="lincdisease", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("fetch", help="All public downloads into --data-dir (needs network); checks hand-placed files")
    s.add_argument("--data-dir", required=True)
    s.set_defaults(func=cmd_fetch)

    s = sub.add_parser("qc", help="Step 1: QC flags, copy-attribution confidence, confidence tier")
    _run_inputs(s)
    s.add_argument("--peptides-dir", required=True,
                   help="directory with proteomics/ and peptidomics/ results of the lincorfz notebook")
    s.set_defaults(func=cmd_qc)

    s = sub.add_parser("homology", help="Step 1: DIAMOND homology to annotated proteins")
    _run_inputs(s)
    s.add_argument("--diamond", default="diamond")
    _threads(s)
    s.set_defaults(func=cmd_homology)

    s = sub.add_parser("regions", help="Step 2: quantifiable exonic regions of loci and controls")
    _run_inputs(s)
    s.add_argument("--min-quant-bases", type=int, default=100)
    s.set_defaults(func=cmd_regions)

    s = sub.add_parser("manifest", help="Step 2: recount3 samples (TCGA-LIHC, GTEx)")
    s.add_argument("--data-dir", required=True)
    # 10 samples per non-liver tissue suffice for the tissue medians of liver and tumour specificity (D6,
    # actionability) and keep streaming time down; the seed fixes which samples are drawn
    s.add_argument("--gtex-per-tissue", type=int, default=10)
    s.add_argument("--seed", type=int, default=1)
    s.set_defaults(func=cmd_manifest)

    s = sub.add_parser("quant", help="Step 2: stream recount3 bigWigs and sum coverage over the regions")
    s.add_argument("--project", nargs="*", help="restrict to these recount3 projects")
    s.add_argument("--worker", type=int, default=0)
    s.add_argument("--n-workers", type=int, default=1)
    s.add_argument("--scratch", help="directory for the temporary bigWig downloads")
    s.add_argument("--sample-set", default="quant",
                   help="subdirectory of --out-dir with samples.tsv (quant: TCGA/GTEx; quant_sra, quant_hcc, "
                        "quant_nmd: SRA cohorts)")
    s.add_argument("--features", nargs="*", help="quantify only these features (added later, e.g. a new control)")
    s.add_argument("--sums-dir", default="sums", help="subdirectory of the sample set for the sums (sums_<name> for "
                                                       "added features; joined to sums/ by sample)")
    s.set_defaults(func=cmd_quant)

    s = sub.add_parser("scaling", help="Step 2: size factors of the TCGA-LIHC / GTEx samples (recount3 gene sums)")
    s.add_argument("--data-dir", required=True)
    s.set_defaults(func=cmd_scaling)

    s = _cohort_step(sub, "chronic", "quant_sra", "Step 3, D1 chronic liver disease (SRA cohorts in quant_sra/)")
    s.add_argument("--govaere-dir", help="Govaere directory (locus_cov/, counts/, star_logs/); its runs enter as "
                                         "the cohort 'govaere'")
    s.add_argument("--govaere-exclude", nargs="*", default=[], help="Govaere runs left out (config govaere.exclude)")
    s.add_argument("--govaere-runs", help="ENA run table of the Govaere runs")
    s.add_argument("--govaere-soft", help="GEO sample annotation of GSE135251 (acc.cgi text format)")
    s.set_defaults(func=cmd_chronic)

    s = _cohort_step(sub, "hcc", "quant_hcc", "Step 3, D1 HCC replication in independent cohorts (quant_hcc/)")
    s.set_defaults(func=cmd_hcc)

    s = sub.add_parser("somatic", help="Step 3, D3: somatic copy number in TCGA-LIHC")
    _run_inputs(s)
    _pseudocount(s)
    s.set_defaults(func=cmd_somatic)

    s = sub.add_parser("genetics", help="Step 3, D4: GWAS Catalog liver traits near each feature, GTEx liver eQTLs")
    _run_inputs(s)
    s.add_argument("--flank", type=int, default=50_000)
    s.set_defaults(func=cmd_genetics)

    s = sub.add_parser("literature", help="Step 3, D5: Europe PMC counts")
    _run_inputs(s)
    _threads(s)
    s.set_defaults(func=cmd_literature)

    s = sub.add_parser("disease", help="Step 3: disease components D1, D2, D6, collects D3-D5, locus disease score")
    _run_inputs(s)
    _pseudocount(s)
    s.set_defaults(func=cmd_disease)

    s = sub.add_parser("weights", help="Fitted disease-score weights (alternative score, reported)")
    _run_inputs(s)
    s.add_argument("--min-papers", type=int, default=weights.MIN_PAPERS,
                   help="liver papers (non-retracted) for an established liver-disease lncRNA")
    s.set_defaults(func=cmd_weights)

    s = sub.add_parser("peptide", help="Step 4 M1 sequence features and actionability HLA-I presentation")
    _run_inputs(s)
    s.add_argument("--tool-env", required=True, help="conda env with hmmer, metapredict and mhcflurry")
    _threads(s)
    s.add_argument("--skip-hla", action="store_true")
    s.set_defaults(func=cmd_peptide)

    s = sub.add_parser("support", help="Step 4 M1: ORF vs exon PhyloP contrast, published CRISPR screen hits")
    _run_inputs(s)
    s.add_argument("--phylop", default=peptide.PHYLOP_URL, help="PhyloP100way bigWig (path or URL)")
    _threads(s)
    s.set_defaults(func=cmd_support)

    s = sub.add_parser("encode", help="Step 4: ENCODE HepG2 fractionation, eCLIP and chromatin")
    _run_inputs(s)
    s.set_defaults(func=cmd_encode)

    s = _cohort_step(sub, "nmd", "quant_nmd", "Step 4, M2: NMD sensitivity from HeLa NMD knock-downs (quant_nmd/)")
    s.set_defaults(func=cmd_nmd)

    s = sub.add_parser("mechanism", help="Step 4: mechanism class features and scores")
    _run_inputs(s)
    _pseudocount(s)
    s.set_defaults(func=cmd_mechanism)

    s = sub.add_parser("report", help="Outputs: per-ORF table, candidates, control panel")
    _run_inputs(s)
    s.set_defaults(func=cmd_report)

    s = sub.add_parser("html", help="Outputs: sortable HTML table of the candidates and screen-supported ORFs")
    s.set_defaults(func=cmd_html)

    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main()
