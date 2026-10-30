"""
Command-line driver of lincdisease; one sub-command per step of LOGIC.md, each writing into --out-dir.

  python -m lincdisease qc       --results-dir R --candidates T --data-dir D --out-dir O
  python -m lincdisease regions  --results-dir R --candidates T --data-dir D --out-dir O
  python -m lincdisease manifest --out-dir O [--gtex-per-tissue 10]
  python -m lincdisease quant    --out-dir O --project LIHC [--worker 0 --n-workers 4] [--scratch S]
  python -m lincdisease literature --results-dir R --candidates T --data-dir D --out-dir O
  python -m lincdisease genetics --data-dir D --out-dir O
  python -m lincdisease mechanism --results-dir R --candidates T --data-dir D --peptides-dir P --out-dir O
  python -m lincdisease disease  --results-dir R --candidates T --data-dir D --out-dir O
  python -m lincdisease report   --results-dir R --candidates T --data-dir D --out-dir O

R: lincorfz results/ directory (background); T: candidate ORF table; D: download cache of public resources.
"""
import argparse
import os

import pandas as pd

import glob

from . import controls, disease, genetics, inputs, integrate, literature, mechanism, qc, quant, report, resources
from . import scoring


def _common(p):
    p.add_argument("--results-dir", required=True, help="lincorfz results/ directory (background)")
    p.add_argument("--candidates", required=True, help="candidate ORF table (lincorfz format)")
    p.add_argument("--data-dir", required=True, help="download cache for public resources")


def cmd_qc(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    kw = dict(segdups=resources.read_segdups(args.data_dir), false_dups=resources.read_false_dups(args.data_dir))
    feats_path = os.path.join(args.out_dir, "quant", "features.tsv")
    if os.path.exists(feats_path):
        feats = pd.read_csv(feats_path, sep="\t", index_col=0)
        kw["quantifiable"] = orfs.set_index("orf_id")["locus_id"].map(feats["quantifiable"])
    res = qc.run_qc(orfs, **kw)
    out = pd.concat([orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate"]],
                     res.drop(columns="orf_id")], axis=1)
    path = os.path.join(args.out_dir, "qc.orfs.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"wrote {path}: {len(out)} ORFs")
    print("tiers, all ORFs:", out["confidence_tier"].value_counts().to_dict())
    print("tiers, candidates:", out.loc[out["is_candidate"], "confidence_tier"].value_counts().to_dict())


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


def cmd_manifest(args):
    """TCGA-LIHC (all samples), GTEx liver (all) and up to --gtex-per-tissue samples of every other GTEx tissue."""
    parts = [quant.manifest("tcga", "LIHC", TCGA_COLUMNS)]
    projects = pd.read_csv(resources.fetch("gtex_projects", args.data_dir), sep="\t")["project"].unique()
    for project in sorted(set(projects) - {"STUDY_NA"}):   # STUDY_NA: tissue unknown
        m = quant.manifest("gtex", project, GTEX_COLUMNS)
        if project != "LIVER" and len(m) > args.gtex_per_tissue:
            m = m.sample(args.gtex_per_tissue, random_state=args.seed)
        parts.append(m)
    man = pd.concat(parts, ignore_index=True)
    d = os.path.join(args.out_dir, "quant")
    os.makedirs(d, exist_ok=True)
    man.to_csv(os.path.join(d, "samples.tsv"), sep="\t", index=False)
    print(man.groupby(["source", "project"]).size().to_string())


def cmd_quant(args):
    d = os.path.join(args.out_dir, "quant")
    regions = pd.read_csv(os.path.join(d, "regions.tsv"), sep="\t")
    samples = pd.read_csv(os.path.join(d, "samples.tsv"), sep="\t", low_memory=False)
    if args.project:
        samples = samples[samples["project"].isin(args.project)]
    samples = samples.iloc[args.worker::args.n_workers]
    path = quant.run(samples, regions, os.path.join(d, "sums"), worker=f"{args.worker}", scratch=args.scratch)
    print(f"done: {path}")


def feature_names(results_dir, candidates):
    """feature_id -> symbol for every locus (its locus_name) and control ('ctrl:<gene>')."""
    loci = inputs.loci(inputs.load(results_dir, candidates))
    names = loci["locus_name"].copy()
    for g in controls.control_genes():
        names[f"ctrl:{g}"] = g
    return names


def cmd_literature(args):
    names = feature_names(args.results_dir, args.candidates)
    aliases = literature.hgnc_aliases(resources.fetch("hgnc", args.data_dir))
    res = literature.d5(names, aliases, os.path.join(args.data_dir, "europepmc_counts.json"), threads=args.threads)
    path = os.path.join(args.out_dir, "D5.literature.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; features with queries failing: {res['n_liver'].isna().sum()}")


def cmd_genetics(args):
    feats = pd.read_csv(os.path.join(args.out_dir, "quant", "features.tsv"), sep="\t", index_col=0)
    assoc = genetics.read_liver_associations(resources.fetch("gwas_catalog", args.data_dir))
    res = genetics.d4(feats[["chrom", "start0", "end1"]], assoc, flank=args.flank)
    path = os.path.join(args.out_dir, "D4.gwas.tsv")
    res.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}; {len(assoc)} liver-trait associations; features with a hit: {(res['gwas_n'] > 0).sum()}")


def cmd_mechanism(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    qct = pd.read_csv(os.path.join(args.out_dir, "qc.orfs.tsv"), sep="\t").set_index("orf_id").reindex(orfs["orf_id"])
    qct.index = orfs.index
    pdir = args.peptides_dir
    read = lambda *parts: pd.read_csv(os.path.join(pdir, *parts), sep="\t")
    peptides = mechanism.specific_peptides(
        read("proteomics", "analysis", "lincORF_peptides.tsv"), read("proteomics", "fasta", "lincORF_entries.tsv"),
        read("peptidomics", "analysis", "immunopeptides.tsv"), read("peptidomics", "fasta", "lincORF_entries.tsv"))
    tx_exons = quant.read_gtf_exons(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.gtf"))
    loci_tx = pd.read_csv(os.path.join(args.results_dir, "lincRNAs", "lincRNA_transcripts.tsv"), sep="\t")
    gencode = quant.read_gtf_exons(resources.fetch("gencode_gtf", args.data_dir))
    feats = mechanism.features(orfs, qct, peptides, mechanism.stop_to_last_junction(orfs, tx_exons),
                               mechanism.neighbour_position(loci_tx, gencode),
                               mechanism.small_rna_hosts(loci_tx, mechanism.small_rna_genes(
                                   gencode, resources.read_wgrna(args.data_dir))))
    res = mechanism.class_table(feats, pd.Series(True, index=orfs.index))
    out = pd.concat([orfs[["orf_id", "orf_label", "locus_id", "locus_name", "is_candidate"]],
                     feats.add_prefix("raw_"), res], axis=1)
    path = os.path.join(args.out_dir, "mechanism.orfs.tsv")
    out.to_csv(path, sep="\t", index=False)
    peptides.to_csv(os.path.join(args.out_dir, "peptides.orfs.tsv"), sep="\t", index_label="orf_id")
    print(f"wrote {path}")
    print("dominant class, candidates:", out.loc[out["is_candidate"], "dominant_class"].value_counts().to_dict())


def load_expression(out_dir, pseudocount):
    d = os.path.join(out_dir, "quant")
    sums = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(d, "sums", "sums.*.tsv")))])
    sums = sums.drop_duplicates("sample_id")
    feats = pd.read_csv(os.path.join(d, "features.tsv"), sep="\t", index_col=0)
    samples = pd.read_csv(os.path.join(d, "samples.tsv"), sep="\t", low_memory=False)
    samples = samples[samples["sample_id"].isin(sums["sample_id"])]
    ncov = quant.normalize(sums, feats, samples)
    return ncov, disease.log_expr(ncov, pseudocount), samples, feats


def cmd_disease(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    loci = inputs.loci(orfs)
    ncov, expr, samples, feats = load_expression(args.out_dir, args.pseudocount)
    print(samples.groupby(["source", "project"]).size().to_string())
    cdr = pd.read_excel(resources.fetch("tcga_cdr", args.data_dir), sheet_name=0)

    d1 = disease.d1_hcc(expr, samples)
    d2 = disease.d2_survival(expr, samples, cdr, ncov=ncov)
    cell = loci[["log2TPM_HepG2", "log2TPM_Huh7"]].astype(float).max(axis=1)
    d6 = integrate.liver_context(expr, samples, cell)
    d4 = pd.read_csv(os.path.join(args.out_dir, "D4.gwas.tsv"), sep="\t", index_col=0)
    d5 = pd.read_csv(os.path.join(args.out_dir, "D5.literature.tsv"), sep="\t", index_col=0)

    stats = pd.DataFrame(index=feats.index)
    stats["D1_HCC"] = d1["D1_HCC"]
    stats["D1_chronic"] = float("nan")           # disease cohorts not yet quantified
    stats["D2"] = d2["D2"]
    stats["D3"] = float("nan")                   # copy number / methylation not yet implemented
    stats["D4"] = d4["D4"]
    stats["D5"] = d5["D5"]
    stats = stats.join(d6)
    # loci that cannot be measured in unstranded public RNA-seq have no expression-based components
    unq = ~feats["quantifiable"].reindex(stats.index).fillna(False).astype(bool)
    stats.loc[unq, ["D1_HCC", "D2", "gtex_liver_median", "liver_specificity"]] = float("nan")
    is_bg = pd.Series(stats.index.isin(loci.index), index=stats.index)
    score = integrate.disease_score(stats, is_bg)

    tum_spec = integrate.tumour_specificity(expr, samples)
    tum_spec[unq.reindex(tum_spec.index).fillna(False)] = float("nan")
    out = pd.concat([stats.add_prefix("stat_"), score, d1.drop(columns="D1_HCC"), d2.drop(columns="D2"),
                     d4.drop(columns="D4"), d5.drop(columns=["D5", "symbol"]), tum_spec.rename("tumour_specificity")], axis=1)
    out.insert(0, "symbol", loci["locus_name"].reindex(out.index).fillna(out.index.to_series().str.replace("ctrl:", "")))
    out.insert(1, "is_control", out.index.str.startswith("ctrl:"))
    out.insert(2, "is_candidate", loci["is_candidate"].reindex(out.index).fillna(False).astype(bool))
    out = out.sort_values("disease_score", ascending=False)
    path = os.path.join(args.out_dir, "disease.loci.tsv")
    out.to_csv(path, sep="\t", index_label="feature_id")
    print(f"wrote {path}")


def cmd_report(args):
    orfs = inputs.load(args.results_dir, args.candidates)
    read = lambda name, **kw: pd.read_csv(os.path.join(args.out_dir, name), sep="\t", **kw)
    qct = read("qc.orfs.tsv").set_index("orf_id").reindex(orfs["orf_id"]).set_axis(orfs.index)
    mech = read("mechanism.orfs.tsv").set_index("orf_id").reindex(orfs["orf_id"]).set_axis(orfs.index)
    dis = read("disease.loci.tsv", index_col=0)
    peptides = read("peptides.orfs.tsv", index_col=0)
    act = report.actionability(orfs, dis, peptides)
    table = report.orf_table(orfs, qct, dis, mech, act)
    table.to_csv(os.path.join(args.out_dir, "lincdisease.orfs.tsv"), sep="\t", index=False)
    cand = table[table["is_candidate"]].sort_values(["disease_score", "priority_score"], ascending=False)
    cand.to_csv(os.path.join(args.out_dir, "lincdisease.candidates.tsv"), sep="\t", index=False)
    by_name = mech.assign(rank=orfs["rank"]).sort_values("rank").drop_duplicates("locus_name").set_index("locus_name")
    ctrl = report.control_report(dis, by_name)
    ctrl.to_csv(os.path.join(args.out_dir, "controls.tsv"), sep="\t", index=False)
    print(ctrl.to_string(index=False))
    print(f"wrote lincdisease.orfs.tsv, lincdisease.candidates.tsv, controls.tsv to {args.out_dir}")


def main():
    p = argparse.ArgumentParser(prog="lincdisease", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("qc", help="Step 1: QC flags, copy-attribution confidence, confidence tier")
    _common(q)
    q.set_defaults(func=cmd_qc)

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
    s.set_defaults(func=cmd_quant)

    lt = sub.add_parser("literature", help="Step 3, D5: Europe PMC counts")
    _common(lt)
    lt.add_argument("--threads", type=int, default=4)
    lt.set_defaults(func=cmd_literature)

    gw = sub.add_parser("genetics", help="Step 3, D4: GWAS Catalog liver traits near each feature")
    gw.add_argument("--data-dir", required=True)
    gw.add_argument("--flank", type=int, default=50_000)
    gw.set_defaults(func=cmd_genetics)

    mc = sub.add_parser("mechanism", help="Step 4: mechanism class features and scores")
    _common(mc)
    mc.add_argument("--peptides-dir", required=True,
                    help="directory with proteomics/ and peptidomics/ results of the lincorfz notebook")
    mc.set_defaults(func=cmd_mechanism)

    ds = sub.add_parser("disease", help="Step 3: disease components D1, D2, D6 and the locus disease score")
    _common(ds)
    ds.add_argument("--pseudocount", type=float, default=0.05, help="added to ncov before log2")
    ds.set_defaults(func=cmd_disease)

    rp = sub.add_parser("report", help="Outputs: per-ORF table, candidates, control panel")
    _common(rp)
    rp.set_defaults(func=cmd_report)

    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    args.func(args)


if __name__ == "__main__":
    main()
