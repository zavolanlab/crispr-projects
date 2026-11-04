#!/usr/bin/env python
"""
Integrate all evidence layers per ORF, rank ORFs and write browser-ready annotations.

Periodicity p-value (per replicate and pooled) = min of two tests:
  * per-codon test of the HCC Ribo-seq paper (one-sided Wilcoxon rank-sum, in-frame vs each off-frame,
    Fisher-combined; computed by call_orfs.py) -- robust, but underpowered for short ORFs whose
    P-sites pile up on a few codons;
  * read-count frame test: one-sided binomial test that the in-frame share of the ORF's P-sites exceeds
    1/3, applied only if frame 0 dominates, >= --frame-test-min-codons codons carry in-frame P-sites and
    the in-frame fraction is >= --frame-test-min-frac-inframe (guards against a single stacked codon).
  Both p-values are reported (<label>_p_periodicity, <label>_p_frame_count).

Bayesian frame test (bayes_frame_test.py, read-length-specific priors): evidence classes
  strong   pooled log10 BF >= --bayes-strong, in-frame P-sites >= --min-inframe-strong and
           >= --min-reps-strong samples with log10 BF >= --bayes-rep-support
  moderate pooled log10 BF >= --bayes-moderate, in-frame P-sites >= --min-inframe-moderate
  weak     pooled log10 BF >= --bayes-weak, in-frame P-sites >= 3
The final class is the higher of the Bayesian class and the periodicity-test class below.

Periodicity-test evidence class (pooled replicates; thresholds are parameters):
  strong    p_periodicity < --p-strong, in-frame P-sites >= --min-inframe-strong,
            in-frame fraction >= --min-frac-inframe-strong, >= --min-reps-strong replicates p < 0.05
  moderate  p < --p-moderate, in-frame P-sites >= --min-inframe-moderate, in-frame fraction >= 0.4
  weak      p < --p-weak, in-frame P-sites >= 3
  none      otherwise

Multimapping reads: the input BAMs have multimappers resolved upstream (one alignment per read,
NH:1, chosen among the tied loci with probabilities from unique-read counts), so all reads count
towards the evidence class. MAPQ keeps the value of the original alignment (< 255 for reads that
had several loci); ORFs with < --min-frac-unique of their P-sites from MAPQ 255 reads are only
flagged in 'resolved_multimapper_dominated' (and in the GTF), not demoted.

Priority score (0..1) = weighted mean of components, each a percentile rank among ORFs with
Ribo-seq evidence (missing values -> 0, except conservation -> 0.5):
  MRL           isoform-level MRL percentile among canonical mRNAs        (--w-mrl)
  ribo          -log10 pooled periodicity p-value                          (--w-ribo)
  island        log2 P-site density ORF / rest of the isoform              (--w-island)
  expression    Sanity log2 TPM of the locus in --primary-condition        (--w-expr)
  conservation  PhyloCSF median                                            (--w-cons)
  uniqueness    fraction of ORF P-sites from MAPQ 255 reads (raw)          (--w-unique, 0 = off)
  peptides      1 if a unique peptide was matched (only if peptides given) (--w-pep)

Outputs (prefix lincORFs):
  .ranked.tsv          all ORFs with every metric, sorted by score
  .loci.ranked.tsv     best ORF per lincRNA locus
  .igv.gtf             expressed lincRNA isoforms + ORFs with >= --gtf-min-evidence (CDS rows)
  .orfs.bed            BED12 of the same ORFs (colour = evidence class, score = priority)
  .faa                 protein sequences of the same ORFs
  .summary.pdf         QC plots of the score components
"""
import argparse
import os
import subprocess
import tempfile

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from lincorfz_utils import read_gtf  # noqa: E402

EVIDENCE_ORDER = {"strong": 4, "moderate": 3, "weak": 2, "none": 0}
EVIDENCE_RGB = {"strong": "178,24,43", "moderate": "239,138,98", "weak": "103,169,207",
                "none": "150,150,150"}
SMALL_RNA_TYPES = {"snoRNA", "snRNA", "miRNA", "misc_RNA", "scaRNA", "rRNA", "rRNA_pseudogene",
                   "Mt_rRNA", "Mt_tRNA", "vault_RNA", "ribozyme", "scRNA", "sRNA", "Y_RNA"}


def orf_blocks(s):
    return [tuple(map(int, b.split("-"))) for b in s.split(";")]


def trim_blocks(blocks, strand, n, end):
    """Split n nt off the 5' ('5p') or 3' ('3p') end of genomic blocks (transcript orientation).
    Returns (terminal_blocks, remaining_blocks); splice-aware."""
    order = list(blocks) if strand == "+" else list(blocks)[::-1]
    if end == "3p":
        order = order[::-1]
    term, rest, left = [], [], n
    for s, e in order:
        ln = e - s
        if left <= 0:
            rest.append((s, e))
            continue
        take = min(left, ln)
        from_low = (strand == "+") == (end == "5p")   # terminal nt at the low-coordinate side?
        if from_low:
            term.append((s, s + take))
            if take < ln:
                rest.append((s + take, e))
        else:
            term.append((e - take, e))
            if take < ln:
                rest.append((s, e - take))
        left -= take
    return sorted(term), sorted(rest)


def bed_overlaps(orfs, ref_df, name_col, strand_flag, tmp, tag):
    a = os.path.join(tmp, f"{tag}_orfs.bed")
    b = os.path.join(tmp, f"{tag}_ref.bed")
    rows = [(r.chrom, s, e, r.orf_id, 0, r.strand) for r in orfs.itertuples() for s, e in orf_blocks(r.orf_blocks)]
    pd.DataFrame(rows).sort_values([0, 1]).to_csv(a, sep="\t", header=False, index=False)
    ref_df[["chrom", "start0", "end1", name_col]].assign(score=0, strand=ref_df["strand"]).sort_values(
        ["chrom", "start0"]).to_csv(b, sep="\t", header=False, index=False)
    out = subprocess.run(["bedtools", "intersect", "-a", a, "-b", b, strand_flag, "-wa", "-wb"],
                         check=True, capture_output=True, text=True).stdout
    res = {}
    for line in out.splitlines():
        f = line.split("\t")
        res.setdefault(f[3], set()).add(f[9])
    return pd.Series({k: ",".join(sorted(v)) for k, v in res.items()}, dtype=object)


def pct_rank(s, fill):
    r = s.rank(pct=True)
    return r.fillna(fill)


def gtf_attr(d):
    return " ".join(f'{k} "{v}";' for k, v in d.items() if v is not None and v == v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orfs", required=True)
    ap.add_argument("--transcripts-tsv", required=True)
    ap.add_argument("--transcripts-gtf", required=True)
    ap.add_argument("--transcript-ribo", required=True)
    ap.add_argument("--polysome-transcripts", required=True)
    ap.add_argument("--polysome-orfs", required=True)
    ap.add_argument("--expression", required=True)
    ap.add_argument("--conservation", required=True)
    ap.add_argument("--peptides", default="")
    ap.add_argument("--bayes", default="", help="orf_frame_bayes.tsv (Bayesian frame test)")
    ap.add_argument("--bayes-strong", type=float, default=2.0, help="min pooled log10 BF for 'strong'")
    ap.add_argument("--bayes-moderate", type=float, default=1.0)
    ap.add_argument("--bayes-weak", type=float, default=0.5)
    ap.add_argument("--bayes-rep-support", type=float, default=0.5,
                    help="per-sample log10 BF counted as replicate support (strong needs --min-reps-strong)")
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--hcc-orfs-gtf", default="", help="HCC Ribo-seq paper nORF GTF (lincRNA ORFs with evidence)")
    ap.add_argument("--replicate-labels", required=True, help="comma-separated P-site replicate labels")
    ap.add_argument("--primary-condition", default="HepG2")
    ap.add_argument("--p-strong", type=float, default=1e-3)
    ap.add_argument("--min-inframe-strong", type=float, default=10)
    ap.add_argument("--min-frac-inframe-strong", type=float, default=0.5)
    ap.add_argument("--min-reps-strong", type=int, default=2)
    ap.add_argument("--p-moderate", type=float, default=1e-2)
    ap.add_argument("--min-inframe-moderate", type=float, default=5)
    ap.add_argument("--p-weak", type=float, default=5e-2)
    ap.add_argument("--min-frac-unique", type=float, default=0.1)
    ap.add_argument("--frame-test-min-codons", type=int, default=3)
    ap.add_argument("--frame-test-min-frac-inframe", type=float, default=0.45)
    ap.add_argument("--frame-test", action="store_true",
                    help="also use the frequentist read-count frame test (superseded by the Bayesian test)")
    ap.add_argument("--w-mrl", type=float, default=0.25)
    ap.add_argument("--w-ribo", type=float, default=0.20)
    ap.add_argument("--w-island", type=float, default=0.10)
    ap.add_argument("--w-expr", type=float, default=0.20)
    ap.add_argument("--w-cons", type=float, default=0.10)
    ap.add_argument("--w-unique", type=float, default=0.15)
    ap.add_argument("--w-pep", type=float, default=0.05)
    # 'multimapping' (former class, now a flag) is accepted as an alias of 'weak'
    ap.add_argument("--gtf-min-evidence", default="weak", choices=list(EVIDENCE_ORDER) + ["multimapping"])
    ap.add_argument("--prefix", default="lincORFs")
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(dir=".")
    orfs = pd.read_csv(args.orfs, sep="\t")
    tx = pd.read_csv(args.transcripts_tsv, sep="\t").set_index("transcript_id")
    reps = [r for r in args.replicate_labels.split(",") if r]

    # ---------------------------------------------------------------- join evidence layers
    ptx = pd.read_csv(args.polysome_transcripts, sep="\t").set_index("transcript_id")
    orfs = orfs.merge(ptx[["MRL", "MRL_pct_vs_mRNA", "polysome_mean_coverage", "frac_polysomal_ge2"]],
                      left_on="transcript_id", right_index=True, how="left")
    porf = pd.read_csv(args.polysome_orfs, sep="\t")
    orfs = orfs.merge(porf[["orf_id", "orf_MRL", "orf_MRL_pct_vs_mRNA", "orf_polysome_mean_coverage"]],
                      on="orf_id", how="left")
    expr = pd.read_csv(args.expression, sep="\t").set_index("gene_id")
    ecols = [c for c in expr.columns if c.startswith(("log2CPM_", "log2TPM_", "err_log2"))
             and not c.startswith("log2CPM_sample_")]
    orfs = orfs.merge(expr[ecols], left_on="locus_id", right_index=True, how="left")
    trib = pd.read_csv(args.transcript_ribo, sep="\t").set_index("transcript_id")
    orfs["tx_psites"] = orfs["transcript_id"].map(trib["tx_psites"])
    cons = pd.read_csv(args.conservation, sep="\t")
    orfs = orfs.merge(cons, on="orf_id", how="left")
    has_bayes = bool(args.bayes) and os.path.exists(args.bayes)
    if has_bayes:
        orfs = orfs.merge(pd.read_csv(args.bayes, sep="\t"), on="orf_id", how="left")
    has_pep = bool(args.peptides) and os.path.exists(args.peptides)
    if has_pep:
        pep = pd.read_csv(args.peptides, sep="\t")
        orfs = orfs.merge(pep, on="orf_id", how="left")
        orfs["n_unique_peptides"] = orfs["n_unique_peptides"].fillna(0).astype(int)
    for c in ["gene_name_source", "name_source_gene_ids", "n_exons", "lncRNA_gene_names", "lncRNA_gene_ids",
              "pseudogene_names", "pseudogene_types",
              "small_RNA_names", "NONCODE_gene_ids", "antisense_coding_exon_overlap",
              "possible_3UTR_extension_of", "possible_5prime_extension_of", "is_reference"]:
        if c in tx.columns and c not in orfs.columns:
            orfs[c] = orfs["transcript_id"].map(tx[c])
    for c in [c for c in tx.columns if c.startswith("TPM_") and c != "TPM_max"]:
        orfs[c] = orfs["transcript_id"].map(tx[c])

    # ORF-level overlaps with reference small RNAs, antisense CDSs and HCC paper ORFs
    ref = read_gtf(args.reference_gtf, features={"exon", "CDS"}, keys=("gene_name", "gene_type"))
    ref = ref.loc[ref["chrom"].isin(set(orfs["chrom"]))]
    small = ref.loc[(ref["feature"] == "exon") & ref["gene_type"].isin(SMALL_RNA_TYPES)]
    orfs["orf_overlaps_small_RNA"] = orfs["orf_id"].map(bed_overlaps(orfs, small, "gene_name", "-s", tmp, "small"))
    cds = ref.loc[ref["feature"] == "CDS"]
    orfs["orf_antisense_to_CDS"] = orfs["orf_id"].map(bed_overlaps(orfs, cds, "gene_name", "-S", tmp, "cds"))
    if args.hcc_orfs_gtf:
        hcc = read_gtf(args.hcc_orfs_gtf, features={"exon"}, keys=("gene_id",))
        hcc = hcc.loc[hcc["gene_id"].str.contains("_lincRNA_1_") & hcc["chrom"].isin(set(orfs["chrom"]))]
        orfs["HCC_paper_lincRNA_ORFs"] = orfs["orf_id"].map(bed_overlaps(orfs, hcc, "gene_id", "-s", tmp, "hcc"))
        # same ORF end (strand-aware) = same frame and stop codon
        ends = {}
        for gid, g in hcc.groupby("gene_id"):
            ends[gid] = (g["chrom"].iat[0], g["strand"].iat[0],
                         g["end1"].max() if g["strand"].iat[0] == "+" else g["start0"].min())
        own_end = [(r.chrom, r.strand, r.orf_end1 if r.strand == "+" else r.orf_start0) for r in orfs.itertuples()]
        orfs["HCC_paper_same_stop"] = [
            any(ends.get(g) == e for g in str(h).split(",")) if isinstance(h, str) else False
            for h, e in zip(orfs["HCC_paper_lincRNA_ORFs"], own_end)]
    subprocess.run(["rm", "-rf", tmp], check=True)

    # ---------------------------------------------------------------- periodicity: Wilcoxon + frame counts
    def frame_count_p(lab):
        f0, f1, f2 = orfs[f"{lab}_f0"], orfs[f"{lab}_f1"], orfs[f"{lab}_f2"]
        n = (f0 + f1 + f2).round().astype(int)
        k = f0.round().astype(int)
        ok = ((f0 > f1) & (f0 > f2) & (orfs[f"{lab}_codons_inframe"] >= args.frame_test_min_codons)
              & (orfs[f"{lab}_frac_inframe"].fillna(0) >= args.frame_test_min_frac_inframe))
        pv = pd.Series(stats.binom.sf(k - 1, n.clip(lower=1), 1.0 / 3.0), index=orfs.index)
        return pv.where(ok & (n > 0), 1.0)

    for lab in reps + ["pooled", "pooled_unique"]:
        orfs[f"{lab}_p_wilcoxon"] = orfs[f"{lab}_p_periodicity"]
        orfs[f"{lab}_p_frame_count"] = frame_count_p(lab)
        if args.frame_test:
            orfs[f"{lab}_p_periodicity"] = np.minimum(orfs[f"{lab}_p_wilcoxon"], orfs[f"{lab}_p_frame_count"])
    orfs["periodicity_test"] = np.where(args.frame_test & (orfs["pooled_p_frame_count"] < orfs["pooled_p_wilcoxon"]),
                                        "frame_count", "wilcoxon")

    # ---------------------------------------------------------------- evidence classes
    p = orfs["pooled_p_periodicity"]
    f0 = orfs["pooled_f0"]
    fr = orfs["pooled_frac_inframe"]
    orfs["n_reps_p05"] = sum((orfs[f"{r}_p_periodicity"] < 0.05).astype(int) for r in reps)
    ev_w = np.where((p < args.p_strong) & (f0 >= args.min_inframe_strong) & (fr >= args.min_frac_inframe_strong)
                    & (orfs["n_reps_p05"] >= args.min_reps_strong), "strong",
                    np.where((p < args.p_moderate) & (f0 >= args.min_inframe_moderate) & (fr >= 0.4), "moderate",
                             np.where((p < args.p_weak) & (f0 >= 3), "weak", "none")))
    orfs["ribo_evidence_periodicity_test"] = ev_w
    order = {"none": 0, "weak": 1, "moderate": 2, "strong": 3}
    if has_bayes:
        bf = orfs["bayes_log10BF_pooled"].fillna(-np.inf)
        orfs["n_reps_bayes_support"] = sum((orfs.get(f"bayes_log10BF_{r}", pd.Series(-np.inf, index=orfs.index))
                                            .fillna(-np.inf) >= args.bayes_rep_support).astype(int) for r in reps)
        ev_b = np.where((bf >= args.bayes_strong) & (f0 >= args.min_inframe_strong)
                        & (orfs["n_reps_bayes_support"] >= args.min_reps_strong), "strong",
                        np.where((bf >= args.bayes_moderate) & (f0 >= args.min_inframe_moderate), "moderate",
                                 np.where((bf >= args.bayes_weak) & (f0 >= 3), "weak", "none")))
        orfs["ribo_evidence_bayes"] = ev_b
        ev = np.where(pd.Series(ev_b).map(order).to_numpy() > pd.Series(ev_w).map(order).to_numpy(), ev_b, ev_w)
    else:
        ev = ev_w
    orfs["ribo_evidence"] = ev
    orfs["resolved_multimapper_dominated"] = (ev != "none") & (orfs["frac_psites_unique"].fillna(0) < args.min_frac_unique)
    uf0 = orfs["pooled_unique_f0"]
    orfs["ribo_evidence_unique_reads"] = np.where(
        (orfs["pooled_unique_p_periodicity"] < args.p_moderate) & (uf0 >= args.min_inframe_moderate), "periodic",
        np.where((orfs["pooled_unique_p_periodicity"] < args.p_weak) & (uf0 >= 3), "weakly_periodic", "none"))

    # ---------------------------------------------------------------- score
    pool = orfs["ribo_evidence"].isin(["strong", "moderate", "weak"])
    prim = f"log2TPM_{args.primary_condition}" if f"log2TPM_{args.primary_condition}" in orfs else f"log2CPM_{args.primary_condition}"
    comp = pd.DataFrame(index=orfs.index)

    def within_pool(series, fill):
        out = pd.Series(np.nan, index=orfs.index)
        out[pool] = pct_rank(series[pool], fill)
        out[~pool] = pct_rank(series, fill)[~pool]
        return out

    comp["c_mrl"] = within_pool(orfs["MRL_pct_vs_mRNA"], 0.0)
    ribo_strength = (orfs["bayes_log10BF_pooled"].clip(lower=-50, upper=50) if has_bayes
                     else -np.log10(p.clip(lower=1e-50)))
    comp["c_ribo"] = within_pool(ribo_strength, 0.0)
    comp["c_island"] = within_pool(orfs["island_log2ratio"], 0.0)
    comp["c_expr"] = within_pool(orfs[prim] if prim in orfs else pd.Series(np.nan, index=orfs.index), 0.0)
    comp["c_cons"] = within_pool(orfs["PhyloCSF_Median"], 0.5)
    comp["c_unique"] = orfs["frac_psites_unique"].fillna(0.0)
    weights = {"c_mrl": args.w_mrl, "c_ribo": args.w_ribo, "c_island": args.w_island, "c_expr": args.w_expr,
               "c_cons": args.w_cons, "c_unique": args.w_unique}
    if has_pep:
        comp["c_pep"] = (orfs["n_unique_peptides"] > 0).astype(float)
        weights["c_pep"] = args.w_pep
    wsum = sum(weights.values())
    orfs = pd.concat([orfs, comp], axis=1)
    orfs["priority_score"] = sum(comp[k] * w for k, w in weights.items()) / wsum
    orfs["evidence_rank"] = orfs["ribo_evidence"].map(EVIDENCE_ORDER)
    orfs = orfs.sort_values(["evidence_rank", "priority_score"], ascending=[False, False]).reset_index(drop=True)
    orfs.insert(0, "rank", np.arange(1, len(orfs) + 1))
    safe = orfs["locus_name"].astype(str).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True).str[:40]
    orfs.insert(1, "orf_label", [f"ORF{r:05d}_{n}_{a}aa" for r, n, a in zip(orfs["rank"], safe, orfs["n_aa"])])

    lead = ["rank", "orf_label", "orf_id", "ribo_evidence", "priority_score", "locus_id", "locus_name",
            "transcript_id", "isoforms", "chrom", "orf_start0", "orf_end1", "strand", "n_aa", "start_codon",
            "kozak_context", "kozak_score", "MRL", "MRL_pct_vs_mRNA", "orf_MRL", prim, "TPM_max",
            "pooled_psites", "pooled_f0", "pooled_frac_inframe", "bayes_log10BF_pooled", "bayes_posterior_pooled",
            "n_reps_bayes_support", "bayes_best_model_pooled", "modal_read_length", "support_read_length_range",
            "support_read_lengths", "read_lengths_with_inframe_evidence", "pooled_p_periodicity",
            "pooled_p_wilcoxon", "pooled_p_frame_count", "pooled_codons_inframe", "n_reps_p05",
            "ribo_evidence_periodicity_test", "ribo_evidence_bayes", "island_log2ratio", "frac_tx_psites_in_orf",
            "frac_psites_unique", "resolved_multimapper_dominated", "ribo_evidence_unique_reads", "gene_name_source",
            "PhyloCSF_Median",
            "PhyloCSF_Power_Median"]
    lead = [c for c in lead if c in orfs.columns]
    orfs = orfs[lead + [c for c in orfs.columns if c not in lead]]
    orfs.to_csv(f"{args.prefix}.ranked.tsv", sep="\t", index=False)
    loci = orfs.drop_duplicates("locus_id")
    loci.to_csv(f"{args.prefix}.loci.ranked.tsv", sep="\t", index=False)
    print(orfs["ribo_evidence"].value_counts().to_string())

    # ---------------------------------------------------------------- GTF / BED / FASTA
    sel = orfs.loc[orfs["evidence_rank"] >= EVIDENCE_ORDER.get(args.gtf_min_evidence, EVIDENCE_ORDER["weak"])]
    tex = read_gtf(args.transcripts_gtf, features={"exon"}, keys=("transcript_id",))
    tex_by = {t: g.sort_values("start0") for t, g in tex.groupby("transcript_id")}
    lines = []

    def add(chrom, s0, e1, feat, strand, attrs, order):
        lines.append((chrom, s0, order, f"{chrom}\tlincorfz\t{feat}\t{s0 + 1}\t{e1}\t.\t{strand}\t.\t{gtf_attr(attrs)}"))

    for tid, g in tex_by.items():
        r = tx.loc[tid]
        a = {"gene_id": r["gene_id"], "transcript_id": tid, "gene_name": r["locus_name"],
             "gene_name_source": r.get("gene_name_source", "NA"),
             "gene_type": "lincRNA_denovo", "transcript_type": "lincRNA_isoform", "TPM_max": f"{r['TPM_max']:.2f}"}
        add(r["chrom"], int(g["start0"].min()), int(g["end1"].max()), "transcript", r["strand"], a, 0)
        for s, e in zip(g["start0"], g["end1"]):
            add(r["chrom"], int(s), int(e), "exon", r["strand"], a, 1)
    for r in sel.itertuples():
        g = tex_by[r.transcript_id]
        a = {"gene_id": r.locus_id, "transcript_id": r.orf_label, "gene_name": r.locus_name,
             "gene_name_source": getattr(r, "gene_name_source", "NA"),
             "gene_type": "lincRNA_denovo", "transcript_type": "lincRNA_ORF", "orf_id": r.orf_id,
             "rank": r.rank, "ribo_evidence": r.ribo_evidence, "priority_score": f"{r.priority_score:.3f}",
             "n_aa": r.n_aa, "start_codon": r.start_codon,
             "bayes_log10BF": f"{r.bayes_log10BF_pooled:.2f}" if has_bayes and pd.notna(r.bayes_log10BF_pooled) else "NA",
             "bayes_posterior": f"{r.bayes_posterior_pooled:.3f}" if has_bayes and pd.notna(r.bayes_posterior_pooled) else "NA",
             "support_read_lengths": getattr(r, "support_read_lengths", "NA") if has_bayes else "NA",
             "p_wilcoxon": f"{r.pooled_p_wilcoxon:.2e}",
             "inframe_psites": int(r.pooled_f0), "codons_inframe": int(r.pooled_codons_inframe),
             "frac_inframe": f"{r.pooled_frac_inframe:.2f}", "island_log2ratio": f"{r.island_log2ratio:.2f}",
             "frac_psites_unique": f"{r.frac_psites_unique:.2f}" if pd.notna(r.frac_psites_unique) else "NA",
             "resolved_multimapper_dominated": str(bool(r.resolved_multimapper_dominated)),
             "MRL": f"{r.MRL:.2f}" if pd.notna(r.MRL) else "NA",
             "MRL_pct_vs_mRNA": f"{r.MRL_pct_vs_mRNA:.2f}" if pd.notna(r.MRL_pct_vs_mRNA) else "NA",
             "expression_" + prim.split("_")[0]: f"{getattr(r, prim):.2f}" if prim in sel and pd.notna(getattr(r, prim)) else "NA",
             "PhyloCSF_median": f"{r.PhyloCSF_Median:.2f}" if pd.notna(r.PhyloCSF_Median) else "NA",
             "aa_seq": r.aa_seq}
        add(r.chrom, int(g["start0"].min()), int(g["end1"].max()), "transcript", r.strand, a, 0)
        for s, e in zip(g["start0"], g["end1"]):
            add(r.chrom, int(s), int(e), "exon", r.strand, a, 1)
        blocks = orf_blocks(r.orf_blocks)
        # GTF convention: CDS excludes the stop codon, which gets its own rows (split if spliced)
        stop_blocks, cds_blocks = trim_blocks(blocks, r.strand, 3, "3p")
        start_blocks, _ = trim_blocks(blocks, r.strand, 3, "5p")
        for s, e in cds_blocks:
            add(r.chrom, s, e, "CDS", r.strand, a, 2)
        for s, e in start_blocks:
            add(r.chrom, s, e, "start_codon", r.strand, a, 3)
        for s, e in stop_blocks:
            add(r.chrom, s, e, "stop_codon", r.strand, a, 3)
    lines.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(f"{args.prefix}.igv.gtf", "w") as fh:
        fh.write("#gtf-version 2.2\n#lincorfz: de novo lincRNA isoforms and Ribo-seq supported ORFs\n")
        fh.write("\n".join(l[3] for l in lines) + "\n")

    with open(f"{args.prefix}.orfs.bed", "w") as fh:
        fh.write('track name="lincORFs" itemRgb="On" description="lincorfz ORFs (colour = Ribo-seq evidence)"\n')
        for r in sel.sort_values(["chrom", "orf_start0"]).itertuples():
            b = orf_blocks(r.orf_blocks)
            fh.write("\t".join(map(str, [r.chrom, b[0][0], b[-1][1], r.orf_label, int(round(1000 * r.priority_score)),
                                         r.strand, b[0][0], b[-1][1], EVIDENCE_RGB[r.ribo_evidence], len(b),
                                         ",".join(str(e - s) for s, e in b) + ",",
                                         ",".join(str(s - b[0][0]) for s, _ in b) + ","])) + "\n")
    with open(f"{args.prefix}.faa", "w") as fh:
        for r in sel.itertuples():
            fh.write(f">{r.orf_label} {r.orf_id} {r.chrom}:{r.orf_start0 + 1}-{r.orf_end1}({r.strand}) "
                     f"evidence={r.ribo_evidence}\n{r.aa_seq}\n")

    # ---------------------------------------------------------------- QC plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, (col, lab) in zip(axes.flat, [("MRL", "isoform MRL"), ("pooled_p_periodicity", "-log10 p periodicity"),
                                          ("island_log2ratio", "island log2 ratio"), (prim, f"log2 CPM {args.primary_condition}"),
                                          ("PhyloCSF_Median", "PhyloCSF median"), ("priority_score", "priority score")]):
        if col not in orfs:
            continue
        for e in ["strong", "moderate", "weak", "none"]:
            v = orfs.loc[orfs["ribo_evidence"] == e, col].dropna()
            if col == "pooled_p_periodicity":
                v = -np.log10(v.clip(lower=1e-50))
            if len(v):
                ax.hist(v, bins=40, density=True, histtype="step", label=f"{e} (n={len(v)})")
        ax.set(xlabel=lab, ylabel="density")
    axes.flat[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.prefix}.summary.pdf")


if __name__ == "__main__":
    main()
