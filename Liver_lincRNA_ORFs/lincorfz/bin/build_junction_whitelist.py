#!/usr/bin/env python
"""
Collect splice junctions from STAR SJ.out.tab files and the reference annotation.

Outputs
  junctions.tsv.gz   all RNA-seq junctions with unique/multi read sums per condition
                     (0-based half-open intron coordinates)
  whitelist.tsv.gz   introns (chrom, start0, end1) that are either annotated in the
                     reference GTF or supported by >= --min-unique-reads unique RNA-seq reads
                     in total; spliced Ribo-seq alignments are kept only if all their
                     introns are in this list.
"""
import argparse

import pandas as pd

from lincorfz_utils import exon_chains, introns_of_chain, read_gtf

STRAND = {0: ".", 1: "+", 2: "-"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rnaseq-tsv", required=True, help="sample sheet with columns sample, condition, sj_tab")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--min-unique-reads", type=int, default=2)
    ap.add_argument("--out-junctions", default="junctions.tsv.gz")
    ap.add_argument("--out-whitelist", default="whitelist.tsv.gz")
    args = ap.parse_args()

    sheet = pd.read_csv(args.rnaseq_tsv, sep="\t")
    parts = []
    for _, row in sheet.iterrows():
        sj = pd.read_csv(row["sj_tab"], sep="\t", header=None,
                         names=["chrom", "first", "last", "strand_code", "motif", "annotated",
                                "unique", "multi", "overhang"])
        sj["start0"] = sj["first"] - 1
        sj["end1"] = sj["last"]
        sj["strand"] = sj["strand_code"].map(STRAND)
        sj["condition"] = row["condition"]
        parts.append(sj[["chrom", "start0", "end1", "strand", "motif", "condition", "unique", "multi"]])
    sj = pd.concat(parts)
    wide = sj.pivot_table(index=["chrom", "start0", "end1", "strand", "motif"], columns="condition",
                          values=["unique", "multi"], aggfunc="sum", fill_value=0)
    wide.columns = [f"{v}_{c}" for v, c in wide.columns]
    wide = wide.reset_index()
    wide["unique_total"] = wide[[c for c in wide.columns if c.startswith("unique_")]].sum(axis=1)
    wide["multi_total"] = wide[[c for c in wide.columns if c.startswith("multi_")]].sum(axis=1)
    wide.to_csv(args.out_junctions, sep="\t", index=False)

    exons = read_gtf(args.gtf, features={"exon"}, keys=("transcript_id",))
    annotated = set()
    for chrom, strand, ex in exon_chains(exons).values():
        for s, e in introns_of_chain(ex):
            annotated.add((chrom, s, e))
    supported = set(map(tuple, wide.loc[wide["unique_total"] >= args.min_unique_reads,
                                        ["chrom", "start0", "end1"]].itertuples(index=False)))
    wl = pd.DataFrame(sorted(annotated | supported), columns=["chrom", "start0", "end1"])
    wl.to_csv(args.out_whitelist, sep="\t", index=False)
    print(f"{len(annotated)} annotated introns, {len(supported)} RNA-seq supported introns, "
          f"{len(wl)} in whitelist")


if __name__ == "__main__":
    main()
