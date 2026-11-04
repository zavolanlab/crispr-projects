#!/usr/bin/env python
"""
Select expressed lincRNA transcripts from a StringTie (merged + re-quantified) assembly.

Every assembled transcript is annotated against the reference (GENCODE) on the SAME strand:
  * overlap with exons of CDS-containing transcripts       -> always excluded (mRNA isoforms,
                                                              UTR extensions, uORF/dORF territory)
  * overlap with exons of any transcript of protein-coding genes (retained introns etc.)
                                                           -> excluded unless --keep-pc-gene-exon-overlap
  * located within a protein-coding gene without exon overlap (same-strand intronic,
    frequently pre-mRNA)                                   -> excluded unless --keep-pc-intronic
  * overlaps with lncRNA, pseudogene and small-RNA genes, and with NONCODE transcripts
                                                           -> reported only (names / ids)
Each selected locus is named after the GENCODE gene with the largest exonic overlap (non-small-RNA
genes preferred), else the NONCODE gene with the largest overlap, else the StringTie gene id
(column gene_name_source).
and additionally: antisense exon overlap with coding genes, and the distance to the nearest
same-strand coding transcript (possible UTR extensions are flagged, not removed).

Expression (StringTie re-estimated TPM, mean per condition) and splice-junction support
(unique STAR junction reads summed over all RNA-seq samples, minimum over the introns of a
transcript) are used for filtering. Within a locus (StringTie gene), isoforms below
--min-isoform-fraction of the locus TPM (in the condition where the isoform is highest) are dropped.

Outputs
  lincRNA_candidates.all.tsv   all non-reference-coding assembled transcripts with annotations
                               and the reason for exclusion (or 'selected')
  lincRNA_transcripts.gtf      selected isoforms (transcript + exon rows), sorted
  lincRNA_transcripts.tsv      selected isoforms with annotations
"""
import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

from lincorfz_utils import introns_of_chain, read_chromosomes, read_gtf

SMALL_RNA_TYPES = {"snoRNA", "snRNA", "miRNA", "misc_RNA", "scaRNA", "rRNA", "rRNA_pseudogene",
                   "Mt_rRNA", "Mt_tRNA", "vault_RNA", "ribozyme", "scRNA", "sRNA", "Y_RNA"}


def bedtools(*args):
    out = subprocess.run(["bedtools", *args], check=True, capture_output=True, text=True).stdout
    return out


def write_bed(df, path, name_col, score_col=None):
    b = pd.DataFrame({"chrom": df["chrom"], "s": df["start0"].astype(int), "e": df["end1"].astype(int),
                      "name": df[name_col], "score": df[score_col] if score_col else 0,
                      "strand": df["strand"]})
    b.sort_values(["chrom", "s", "e"]).to_csv(path, sep="\t", header=False, index=False)


def overlap_names(a_bed, b_bed, strand_flag="-s"):
    """dict a_name -> sorted unique b_names overlapping (by exon)."""
    out = bedtools("intersect", "-a", a_bed, "-b", b_bed, strand_flag, "-wa", "-wb")
    res = {}
    for line in out.splitlines():
        f = line.split("\t")
        res.setdefault(f[3], set()).add(f[9])
    return {k: ",".join(sorted(v)) for k, v in res.items()}


def overlap_bp(a_bed, b_bed, strand_flag="-s"):
    """DataFrame (a_name, b_name, bp): exonic overlap in nt of every A x B pair."""
    out = bedtools("intersect", "-a", a_bed, "-b", b_bed, strand_flag, "-wo")
    rows = [(f[3], f[9], int(f[-1])) for f in (l.split("\t") for l in out.splitlines())]
    return pd.DataFrame(rows, columns=["transcript_id", "name", "bp"])


def locus_names(sel, ov_gencode, ov_noncode):
    """
    gene_name per locus (StringTie gene of the selected isoforms):
      1. GENCODE gene with the largest exonic same-strand overlap over the locus' isoforms
         (non-small-RNA genes first; small-RNA genes only if nothing else overlaps),
      2. otherwise the NONCODE gene with the largest overlap,
      3. otherwise the StringTie gene id.
    """
    t2l = sel.set_index("transcript_id")["gene_id"]
    res = {}
    for src, ov in [("GENCODE", ov_gencode), ("NONCODE", ov_noncode)]:
        if ov is None or ov.empty:
            continue
        o = ov.loc[ov["transcript_id"].isin(t2l.index)].copy()
        o["locus"] = o["transcript_id"].map(t2l)
        if src == "GENCODE":
            o["small"] = o["name"].str.split("|").str[2].isin(SMALL_RNA_TYPES)
            o = o.sort_values(["locus", "small", "bp"], ascending=[True, True, False])
            agg = o.groupby(["locus", "name", "small"], as_index=False)["bp"].sum().sort_values(
                ["locus", "small", "bp"], ascending=[True, True, False])
        else:
            agg = o.groupby(["locus", "name"], as_index=False)["bp"].sum().sort_values(
                ["locus", "bp"], ascending=[True, False])
        for locus, g in agg.groupby("locus"):
            if locus in res:
                continue
            top = g["name"].iat[0]
            res[locus] = (top.split("|")[0], src, ",".join(dict.fromkeys(n.split("|")[1] for n in g["name"])),
                          ",".join(dict.fromkeys(n.split("|")[0] for n in g["name"])))
    out = pd.DataFrame.from_dict(res, orient="index", columns=["locus_name", "gene_name_source",
                                                              "name_source_gene_ids", "name_source_gene_names"])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-gtf", required=True)
    ap.add_argument("--quant-gtfs", nargs="+", required=True, help="per-sample stringtie -e GTFs")
    ap.add_argument("--rnaseq-tsv", required=True, help="sample sheet (sample, condition, ...)")
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--junctions", required=True, help="junctions.tsv.gz from build_junction_whitelist.py")
    ap.add_argument("--noncode-gtf", default="")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--chromosomes", default="all")
    ap.add_argument("--min-tx-length", type=int, default=200)
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--min-tpm-monoexonic", type=float, default=5.0)
    ap.add_argument("--min-junction-reads", type=int, default=3)
    ap.add_argument("--min-isoform-fraction", type=float, default=0.1)
    ap.add_argument("--utr-extension-distance", type=int, default=1000)
    ap.add_argument("--keep-pc-gene-exon-overlap", action="store_true")
    ap.add_argument("--keep-pc-intronic", action="store_true")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tmp = tempfile.mkdtemp(dir=args.outdir)
    chroms = set(read_chromosomes(args.chromosomes, args.chrom_sizes)) - {"chrM"}

    # ------------------------------------------------------------------ assembly
    asm = read_gtf(args.merged_gtf, features={"transcript", "exon"})
    asm = asm.loc[asm["chrom"].isin(chroms) & asm["strand"].isin(["+", "-"])]
    tx = asm.loc[asm["feature"] == "transcript"].copy()
    ex = asm.loc[asm["feature"] == "exon"].copy()
    ex["len"] = ex["end1"] - ex["start0"]
    agg = ex.groupby("transcript_id").agg(n_exons=("len", "size"), tx_length=("len", "sum"))
    tx = tx.merge(agg, left_on="transcript_id", right_index=True)

    # expression: TPM per sample from stringtie -e outputs, mean per condition
    sheet = pd.read_csv(args.rnaseq_tsv, sep="\t")
    tpm = {}
    for path in args.quant_gtfs:
        sample = os.path.basename(path).split(".")[0]
        q = read_gtf(path, features={"transcript"}, keys=("transcript_id", "TPM"))
        tpm[sample] = q.set_index("transcript_id")["TPM"].astype(float)
    tpm = pd.DataFrame(tpm).fillna(0.0)
    conds = sheet.drop_duplicates("sample").set_index("sample")["condition"]
    cond_tpm = pd.DataFrame({c: tpm[[s for s in conds.index[conds == c] if s in tpm.columns]].mean(axis=1)
                             for c in conds.unique()})
    cond_tpm.columns = [f"TPM_{c}" for c in cond_tpm.columns]
    tx = tx.merge(cond_tpm, left_on="transcript_id", right_index=True, how="left")
    tpm_cols = list(cond_tpm.columns)
    tx[tpm_cols] = tx[tpm_cols].fillna(0.0)
    tx["TPM_max"] = tx[tpm_cols].max(axis=1)
    tx["best_condition"] = tx[tpm_cols].idxmax(axis=1).str.replace("TPM_", "", regex=False)

    # isoform fraction within StringTie gene (locus), in the isoform's best condition
    for c in tpm_cols:
        tx[f"locus_{c}"] = tx.groupby("gene_id")[c].transform("sum")
    tx["isoform_fraction"] = [
        (r[f"TPM_{r['best_condition']}"] / r[f"locus_TPM_{r['best_condition']}"])
        if r[f"locus_TPM_{r['best_condition']}"] > 0 else 0.0 for _, r in tx.iterrows()]

    # ------------------------------------------------------------------ reference
    ref = read_gtf(args.reference_gtf, features={"exon", "CDS", "gene"},
                   keys=("gene_id", "transcript_id", "gene_name", "gene_type", "transcript_type"))
    ref = ref.loc[ref["chrom"].isin(chroms)]
    cds_tx = set(ref.loc[ref["feature"] == "CDS", "transcript_id"])
    rex = ref.loc[ref["feature"] == "exon"].copy()
    rgene = ref.loc[ref["feature"] == "gene"].copy()

    a_bed = os.path.join(tmp, "asm_exons.bed")
    write_bed(ex.loc[ex["transcript_id"].isin(set(tx["transcript_id"]))], a_bed, "transcript_id")
    a_span = os.path.join(tmp, "asm_spans.bed")
    write_bed(tx, a_span, "transcript_id")

    beds = {}

    def ref_bed(name, df, col="gene_name"):
        p = os.path.join(tmp, f"{name}.bed")
        write_bed(df, p, col)
        beds[name] = p

    ref_bed("coding_exons", rex.loc[rex["transcript_id"].isin(cds_tx)])
    ref_bed("pc_gene_exons", rex.loc[rex["gene_type"] == "protein_coding"])
    ref_bed("pc_gene_spans", rgene.loc[rgene["gene_type"] == "protein_coding"])
    ref_bed("lnc_exons", rex.loc[rex["gene_type"] == "lncRNA"])
    ref_bed("lnc_gene_ids", rex.loc[rex["gene_type"] == "lncRNA"], col="gene_id")
    ref_bed("pseudo_exons", rex.loc[rex["gene_type"].str.contains("pseudogene", na=False)
                                    & ~rex["gene_type"].isin(SMALL_RNA_TYPES)])
    ref_bed("pseudo_types", rex.loc[rex["gene_type"].str.contains("pseudogene", na=False)
                                    & ~rex["gene_type"].isin(SMALL_RNA_TYPES)], col="gene_type")
    ref_bed("small_exons", rex.loc[rex["gene_type"].isin(SMALL_RNA_TYPES)])
    named = rex.loc[rex["gene_type"] != "protein_coding"].copy()
    named["composite"] = named["gene_name"].fillna(named["gene_id"]) + "|" + named["gene_id"] + "|" + named["gene_type"]
    ref_bed("named_exons", named, col="composite")
    cds_spans = rex.loc[rex["transcript_id"].isin(cds_tx)].groupby("transcript_id").agg(
        chrom=("chrom", "first"), start0=("start0", "min"), end1=("end1", "max"),
        strand=("strand", "first"), gene_name=("gene_name", "first")).reset_index()
    ref_bed("coding_tx_spans", cds_spans)

    ann = pd.DataFrame(index=tx["transcript_id"])
    ann["coding_exon_overlap"] = pd.Series(overlap_names(a_bed, beds["coding_exons"]))
    ann["pc_gene_exon_overlap"] = pd.Series(overlap_names(a_bed, beds["pc_gene_exons"]))
    ann["pc_gene_span_overlap"] = pd.Series(overlap_names(a_span, beds["pc_gene_spans"]))
    ann["lncRNA_gene_names"] = pd.Series(overlap_names(a_bed, beds["lnc_exons"]))
    ann["lncRNA_gene_ids"] = pd.Series(overlap_names(a_bed, beds["lnc_gene_ids"]))
    ann["pseudogene_names"] = pd.Series(overlap_names(a_bed, beds["pseudo_exons"]))
    ann["pseudogene_types"] = pd.Series(overlap_names(a_bed, beds["pseudo_types"]))
    ann["small_RNA_names"] = pd.Series(overlap_names(a_bed, beds["small_exons"]))
    ann["antisense_coding_exon_overlap"] = pd.Series(overlap_names(a_bed, beds["coding_exons"], "-S"))
    if args.noncode_gtf:
        nc = read_gtf(args.noncode_gtf, features={"exon"}, keys=("gene_id", "transcript_id"))
        nc = nc.loc[nc["chrom"].isin(chroms) & nc["strand"].isin(["+", "-"])]
        ref_bed("noncode_exons", nc, col="gene_id")
        ann["NONCODE_gene_ids"] = pd.Series(overlap_names(a_bed, beds["noncode_exons"]))
        ov_noncode = overlap_bp(a_bed, beds["noncode_exons"])
        ov_noncode["name"] = ov_noncode["name"] + "|" + ov_noncode["name"]
    else:
        ann["NONCODE_gene_ids"] = np.nan
        ov_noncode = None
    ov_gencode = overlap_bp(a_bed, beds["named_exons"])

    # nearest same-strand coding transcript (signed w.r.t. the lincRNA strand: <0 upstream)
    out = bedtools("closest", "-a", a_span, "-b", beds["coding_tx_spans"], "-s", "-D", "a", "-t", "first")
    near = {}
    for line in out.splitlines():
        f = line.split("\t")
        if f[6] != ".":
            near[f[3]] = (f[9], int(f[-1]))
    ann["nearest_coding_gene_same_strand"] = pd.Series({k: v[0] for k, v in near.items()})
    ann["nearest_coding_distance"] = pd.Series({k: v[1] for k, v in near.items()})
    tx = tx.merge(ann, left_on="transcript_id", right_index=True, how="left")

    # ------------------------------------------------------------------ junction support
    jn = pd.read_csv(args.junctions, sep="\t")
    jkey = dict(zip(zip(jn["chrom"], jn["start0"], jn["end1"]), jn["unique_total"]))
    exg = ex.sort_values("start0").groupby("transcript_id")
    min_sup, n_intr = {}, {}
    for tid, g in exg:
        intr = introns_of_chain(list(zip(g["start0"], g["end1"])))
        c = g["chrom"].iat[0]
        sup = [jkey.get((c, s, e), 0) for s, e in intr]
        min_sup[tid] = min(sup) if sup else np.nan
        n_intr[tid] = len(intr)
    tx["min_intron_junction_reads"] = tx["transcript_id"].map(min_sup)

    # ------------------------------------------------------------------ flags & filters
    d = tx["nearest_coding_distance"]
    tx["possible_3UTR_extension_of"] = np.where(
        (d < 0) & (d.abs() <= args.utr_extension_distance), tx["nearest_coding_gene_same_strand"], "")
    tx["possible_5prime_extension_of"] = np.where(
        (d > 0) & (d <= args.utr_extension_distance), tx["nearest_coding_gene_same_strand"], "")
    tx["pc_intronic_same_strand"] = tx["pc_gene_span_overlap"].notna() & tx["pc_gene_exon_overlap"].isna()
    tx["is_reference"] = ~tx["transcript_id"].str.startswith("MSTRG")

    reason = pd.Series("selected", index=tx.index)

    def flag(mask, why):
        reason.loc[mask & (reason == "selected")] = why

    flag(tx["transcript_id"].isin(cds_tx) | tx["coding_exon_overlap"].notna(), "coding_exon_overlap")
    if not args.keep_pc_gene_exon_overlap:
        flag(tx["pc_gene_exon_overlap"].notna(), "protein_coding_gene_exon_overlap")
    if not args.keep_pc_intronic:
        flag(tx["pc_intronic_same_strand"], "protein_coding_intronic_same_strand")
    flag(tx["tx_length"] < args.min_tx_length, "short")
    flag(tx["TPM_max"] < args.min_tpm, "low_expression")
    flag((tx["n_exons"] == 1) & (tx["TPM_max"] < args.min_tpm_monoexonic), "low_expression_monoexonic")
    flag((tx["n_exons"] > 1) & (tx["min_intron_junction_reads"] < args.min_junction_reads),
         "weak_junction_support")
    flag(tx["isoform_fraction"] < args.min_isoform_fraction, "minor_isoform")
    tx["selection"] = reason

    # names of the selected loci (GENCODE > NONCODE > StringTie id); other transcripts get their
    # own transcript-level GENCODE/NONCODE names for the candidate report
    selected = tx.loc[tx["selection"] == "selected"]
    ln = locus_names(selected, ov_gencode, ov_noncode)
    sel_mask = tx["selection"] == "selected"
    tx["locus_name"] = tx["gene_id"]
    tx["gene_name_source"] = "StringTie"
    tx["name_source_gene_ids"] = ""
    tx["name_source_gene_names"] = ""
    for c in ["locus_name", "gene_name_source", "name_source_gene_ids", "name_source_gene_names"]:
        tx.loc[sel_mask, c] = tx.loc[sel_mask, "gene_id"].map(ln[c]).fillna(tx.loc[sel_mask, c])
    tx.loc[sel_mask & tx["gene_name_source"].eq("StringTie"), "locus_name"] = tx["gene_id"]

    keep_cols = ["gene_id", "transcript_id", "locus_name", "gene_name_source", "name_source_gene_ids",
                 "name_source_gene_names", "chrom", "start0", "end1", "strand", "n_exons",
                 "tx_length", *tpm_cols, "TPM_max", "best_condition", "isoform_fraction",
                 "min_intron_junction_reads", "is_reference", "ref_gene_id", "ref_gene_name",
                 "lncRNA_gene_names", "lncRNA_gene_ids", "pseudogene_names", "pseudogene_types", "small_RNA_names",
                 "NONCODE_gene_ids", "antisense_coding_exon_overlap", "nearest_coding_gene_same_strand",
                 "nearest_coding_distance", "possible_3UTR_extension_of", "possible_5prime_extension_of",
                 "pc_intronic_same_strand", "coding_exon_overlap", "pc_gene_exon_overlap", "selection"]
    keep_cols = [c for c in keep_cols if c in tx.columns]
    report = tx.loc[~tx["transcript_id"].isin(cds_tx), keep_cols].sort_values(["chrom", "start0"])
    report.to_csv(os.path.join(args.outdir, "lincRNA_candidates.all.tsv"), sep="\t", index=False)
    sel = report.loc[report["selection"] == "selected"].copy()
    sel.to_csv(os.path.join(args.outdir, "lincRNA_transcripts.tsv"), sep="\t", index=False)
    print(report["selection"].value_counts().to_string())

    # ------------------------------------------------------------------ GTF of selected isoforms
    sel_ids = set(sel["transcript_id"])
    info = sel.set_index("transcript_id")
    lines = []
    for tid, g in ex.loc[ex["transcript_id"].isin(sel_ids)].groupby("transcript_id"):
        r = info.loc[tid]
        attrs = (f'gene_id "{r["gene_id"]}"; transcript_id "{tid}"; gene_name "{r["locus_name"]}"; '
                 f'gene_name_source "{r["gene_name_source"]}"; gene_type "lincRNA_denovo"; '
                 f'TPM_max "{r["TPM_max"]:.2f}";')
        g = g.sort_values("start0")
        lines.append((r["chrom"], int(g["start0"].min()), 0,
                      f'{r["chrom"]}\tlincorfz\ttranscript\t{g["start0"].min() + 1}\t{g["end1"].max()}\t.\t'
                      f'{r["strand"]}\t.\t{attrs}'))
        exons = list(zip(g["start0"], g["end1"]))
        order = range(len(exons)) if r["strand"] == "+" else range(len(exons) - 1, -1, -1)
        for n, i in enumerate(order, 1):
            s, e = exons[i]
            lines.append((r["chrom"], int(s), 1, f'{r["chrom"]}\tlincorfz\texon\t{s + 1}\t{e}\t.\t{r["strand"]}\t.\t'
                                                  f'{attrs} exon_number "{n}";'))
    lines.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(os.path.join(args.outdir, "lincRNA_transcripts.gtf"), "w") as fh:
        fh.write("\n".join(l[3] for l in lines) + "\n")
    subprocess.run(["rm", "-rf", tmp], check=True)


if __name__ == "__main__":
    main()
