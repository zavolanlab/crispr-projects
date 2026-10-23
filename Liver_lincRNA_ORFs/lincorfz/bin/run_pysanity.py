#!/usr/bin/env python
"""
Normalize RNA-seq gene counts with pySanity (zavolab_pyutils) and report lincRNA loci.

The count matrix combines reference genes (featureCounts on the reference GTF) with the
de novo lincRNA loci (featureCounts on the selected lincRNA isoforms, gene_id = locus),
so that Sanity's library-size model is estimated on the full transcriptome.

Outputs
  sanity.log2cpm.per_sample.tsv   log2 CPM-like expression per sample (all genes)
  sanity.condition_means.tsv      condition means, relative and absolute errors (all genes)
  lincRNA_loci.expression.tsv     the same, restricted to lincRNA loci, plus raw counts
"""
import argparse

import numpy as np
import pandas as pd
from zavolab_pyutils import apply_sanity_normalization_full_bayesian


def read_featurecounts(path):
    fc = pd.read_csv(path, sep="\t", comment="#")
    fc = fc.set_index("Geneid").drop(columns=["Chr", "Start", "End", "Strand", "Length"])
    # column names are BAM paths; the sample sheet maps them to sample ids
    return fc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference-counts", required=True)
    ap.add_argument("--lincrna-counts", required=True)
    ap.add_argument("--rnaseq-tsv", required=True)
    ap.add_argument("--min-total-counts", type=int, default=10,
                    help="reference genes with fewer counts over all samples are dropped")
    ap.add_argument("--n-cores", type=int, default=4)
    args = ap.parse_args()

    sheet = pd.read_csv(args.rnaseq_tsv, sep="\t").drop_duplicates("sample")
    # featureCounts columns are the (staged) BAM names: either the original basename
    # or '<sample>.bam' for chromosome-subset BAMs
    bam2sample = {}
    for _, r in sheet.iterrows():
        bam2sample[r["bam"].split("/")[-1]] = r["sample"]
        bam2sample[f"{r['sample']}.bam"] = r["sample"]

    def rename(df):
        return df.rename(columns=lambda c: bam2sample.get(c.split("/")[-1], c))

    ref = rename(read_featurecounts(args.reference_counts))
    linc = rename(read_featurecounts(args.lincrna_counts))
    samples = list(sheet["sample"])
    ref, linc = ref[samples], linc[samples]
    ref = ref.loc[ref.sum(axis=1) >= args.min_total_counts]
    linc.index = "LOCUS:" + linc.index.astype(str)
    counts = pd.concat([ref, linc]).astype(int)
    print(f"{len(ref)} reference genes, {len(linc)} lincRNA loci, {len(samples)} samples")

    meta = sheet[["sample", "condition"]]
    (per_sample, means, err_rel, err_abs, vg, median_lib, _) = apply_sanity_normalization_full_bayesian(
        counts, meta, sample_col="sample", cond_col="condition", n_cores=args.n_cores)

    per_sample.to_csv("sanity.log2cpm.per_sample.tsv", sep="\t", index_label="gene_id")
    cm = pd.concat([means.add_prefix("log2CPM_"), err_rel.add_prefix("err_rel_"),
                    err_abs.add_prefix("err_abs_")], axis=1)
    cm.to_csv("sanity.condition_means.tsv", sep="\t", index_label="gene_id")

    li = [g for g in counts.index if g.startswith("LOCUS:")]
    out = pd.concat([cm.reindex(li), per_sample.reindex(li).add_prefix("log2CPM_sample_"),
                     counts.loc[li].add_prefix("counts_")], axis=1)
    out.index = [g.replace("LOCUS:", "", 1) for g in out.index]
    # genes without counts are absent from Sanity output -> NaN expression
    out.to_csv("lincRNA_loci.expression.tsv", sep="\t", index_label="gene_id")
    print(out.describe().T.to_string())


if __name__ == "__main__":
    main()
