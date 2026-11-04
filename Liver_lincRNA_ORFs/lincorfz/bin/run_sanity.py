#!/usr/bin/env python
"""
Normalize RNA-seq gene counts with Sanity (Sanity_utils) and report lincRNA loci.

The count matrix combines reference genes (featureCounts on the reference GTF) with the
de novo lincRNA loci (featureCounts on the selected lincRNA isoforms, gene_id = locus), so that
Sanity's normalization is estimated on the full transcriptome. Sanity is run twice:
  * log2 CPM-like expression (no length correction);
  * log2 TPM-like expression (length-corrected with the featureCounts exonic lengths).
Condition means and errors come from ISMARA replicate averaging (Sanity_utils.average_conditions).

Outputs
  sanity.log2cpm.per_sample.tsv, sanity.log2tpm.per_sample.tsv   all genes, per sample
  sanity.condition_means.tsv                                    all genes, per condition
  lincRNA_loci.expression.tsv                                   lincRNA loci: condition means and
                                                                errors (CPM, TPM), per-sample values, counts
"""
import argparse

import pandas as pd
import sanity_utils as su


def read_featurecounts(path):
    fc = pd.read_csv(path, sep="\t", comment="#")
    lengths = fc.set_index("Geneid")["Length"]
    fc = fc.set_index("Geneid").drop(columns=["Chr", "Start", "End", "Strand", "Length"])
    return fc, lengths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference-counts", required=True)
    ap.add_argument("--lincrna-counts", required=True)
    ap.add_argument("--rnaseq-tsv", required=True)
    ap.add_argument("--min-total-counts", type=int, default=10,
                    help="reference genes with fewer counts over all samples are dropped")
    ap.add_argument("--v-method", default="MAP", choices=["MAP", "EAP", "MLE", "MARG"])
    ap.add_argument("--n-threads", type=int, default=4)
    args = ap.parse_args()

    sheet = pd.read_csv(args.rnaseq_tsv, sep="\t").drop_duplicates("sample")
    # featureCounts columns are the (staged) BAM names: the original basename or '<sample>.bam'
    bam2sample = {}
    for _, r in sheet.iterrows():
        bam2sample[r["bam"].split("/")[-1]] = r["sample"]
        bam2sample[f"{r['sample']}.bam"] = r["sample"]

    def rename(df):
        return df.rename(columns=lambda c: bam2sample.get(c.split("/")[-1], c))

    ref, ref_len = read_featurecounts(args.reference_counts)
    linc, linc_len = read_featurecounts(args.lincrna_counts)
    samples = list(sheet["sample"])
    ref, linc = rename(ref)[samples], rename(linc)[samples]
    ref = ref.loc[ref.sum(axis=1) >= args.min_total_counts]
    linc.index = "LOCUS:" + linc.index.astype(str)
    linc_len.index = "LOCUS:" + linc_len.index.astype(str)
    counts = pd.concat([ref, linc]).astype(int)
    lengths = pd.concat([ref_len, linc_len]).reindex(counts.index)
    print(f"{len(ref)} reference genes, {len(linc)} lincRNA loci, {len(samples)} samples")

    meta = sheet[["sample", "condition"]]
    res_cpm = su.sanity_normalization(counts, v_method=args.v_method, n_threads=args.n_threads)
    res_tpm = su.sanity_normalization(counts, v_method=args.v_method, n_threads=args.n_threads,
                                      gene_lengths=lengths)
    avg_cpm = res_cpm.average_conditions(meta)
    avg_tpm = res_tpm.average_conditions(meta)

    res_cpm.log2_expression.to_csv("sanity.log2cpm.per_sample.tsv", sep="\t", index_label="gene_id")
    res_tpm.log2_expression.to_csv("sanity.log2tpm.per_sample.tsv", sep="\t", index_label="gene_id")
    cm = pd.concat([avg_cpm.means.add_prefix("log2CPM_"), avg_cpm.errors.add_prefix("err_log2CPM_"),
                    avg_tpm.means.add_prefix("log2TPM_"), avg_tpm.errors.add_prefix("err_log2TPM_"),
                    res_cpm.gene_mean.rename("log2CPM_mean_all_samples")], axis=1)
    cm.to_csv("sanity.condition_means.tsv", sep="\t", index_label="gene_id")

    li = [g for g in counts.index if g.startswith("LOCUS:")]
    out = pd.concat([cm.reindex(li), res_cpm.log2_expression.reindex(li).add_prefix("log2CPM_sample_"),
                     counts.loc[li].add_prefix("counts_"), lengths.reindex(li).rename("exonic_length")], axis=1)
    out.index = [g.replace("LOCUS:", "", 1) for g in out.index]
    # loci without counts are dropped by Sanity -> NaN expression
    out.to_csv("lincRNA_loci.expression.tsv", sep="\t", index_label="gene_id")
    print(out.describe().T.to_string())
    print("Sanity:", res_cpm.params.get("sanity_version", ""), res_cpm.params.get("sanity_commit", ""))


if __name__ == "__main__":
    main()
