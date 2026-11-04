"""
Peptide homology (LOGIC.md, Step 1).

Homology. Every ORF peptide is searched with DIAMOND blastp (--more-sensitive, default low-complexity masking)
against the GENCODE v42 protein-coding translations and the UniProt human reference proteome (canonical and
additional isoforms). Hits to the ORF's own locus (a GENCODE gene, or the gene named by a UniProt entry,
whose span overlaps the ORF) are discarded. A homolog is the best remaining hit (bitscore) with
E <= `max_evalue` and query coverage >= `min_qcov` %; reported as parent gene, protein id, identity (%),
query coverage (%), parent coverage (%) and the number of distinct homologous genes.

ORF integrity since duplication = parent coverage of the best hit (%): the fraction of the parent protein the
ORF still encodes before its stop codon / frameshift.
"""
import gzip
import os
import re
import subprocess

import numpy as np
import pandas as pd

DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send",
                  "evalue", "bitscore", "qcovhsp", "scovhsp", "qseq_gapped", "sseq_gapped"]


# ---------------------------------------------------------------- sequences

def read_fasta(path):
    """Dict header (without '>') -> sequence; gzip-aware."""
    opener = gzip.open if path.endswith(".gz") else open
    out, name, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def reference_proteins(gencode_translations, uniprot_fastas):
    """
    Reference protein FASTA records for DIAMOND: ids 'gc|<ENST>|<gene>' (GENCODE) and 'up|<accession>|<gene>'
    (UniProt; gene from GN=, else the accession). Returns dict id -> sequence.
    """
    ref = {}
    for h, s in read_fasta(gencode_translations).items():
        f = h.split("|")
        ref[f"gc|{f[1]}|{f[6]}"] = s
    for path in uniprot_fastas:
        for h, s in read_fasta(path).items():
            acc = h.split("|")[1]
            gn = re.search(r" GN=(\S+)", h)
            ref[f"up|{acc}|{gn.group(1) if gn else acc}"] = s
    return ref


def write_fasta(records, path):
    with open(path, "w") as fh:
        for k, s in records.items():
            fh.write(f">{k}\n{s}\n")


# ---------------------------------------------------------------- search

def run_diamond(query_fa, ref_fa, out_tsv, workdir, diamond="diamond", threads=4, max_evalue=1e-3):
    """DIAMOND blastp (--more-sensitive) of query vs reference; returns the hit table."""
    db = os.path.join(workdir, "reference")
    subprocess.run([diamond, "makedb", "--in", ref_fa, "-d", db, "--quiet", "--threads", str(threads)], check=True)
    subprocess.run([diamond, "blastp", "-q", query_fa, "-d", db, "-o", out_tsv, "--more-sensitive", "--quiet",
                    "--threads", str(threads), "--evalue", str(max_evalue), "--max-target-seqs", "50",
                    "--outfmt", "6", *DIAMOND_FIELDS], check=True)
    return read_hits(out_tsv)


def read_hits(path):
    h = pd.read_csv(path, sep="\t", header=None, names=DIAMOND_FIELDS)
    sid = h["sseqid"].str.split("|", expand=True)
    h["source"] = sid[0].map({"gc": "gencode", "up": "uniprot"})
    h["parent_id"] = sid[1]
    h["parent"] = sid[2]
    return h


def gene_spans(gencode_exons):
    """GENCODE gene spans by gene_name: chrom/start0/end1 (all genes of that name)."""
    return gencode_exons.groupby(["gene_name", "chrom"]).agg(start0=("start0", "min"),
                                                            end1=("end1", "max")).reset_index()


def drop_self_hits(hits, orfs, spans):
    """Remove hits whose parent gene overlaps the ORF span (the ORF's own locus); returns (kept, n_self)."""
    # _row identifies a hit: the merge with the gene spans can repeat it
    hits = hits.reset_index(drop=True).assign(_row=lambda d: np.arange(len(d)))
    o = orfs.set_index("orf_id")[["chrom", "orf_start0", "orf_end1"]]
    m = hits.join(o, on="qseqid").merge(spans, left_on=["parent", "chrom"], right_on=["gene_name", "chrom"],
                                        how="left")
    is_self = (m["start0"] < m["orf_end1"]) & (m["end1"] > m["orf_start0"])
    self_rows = set(m.loc[is_self, "_row"])
    keep = ~hits["_row"].isin(self_rows)
    n_self = hits.loc[~keep].groupby("qseqid").size()
    return hits[keep].drop(columns="_row"), n_self


# ---------------------------------------------------------------- per-ORF table

def homology_table(orfs, hits, spans, max_evalue=1e-5, min_qcov=30.0):
    """
    One row per ORF: best non-self homolog (parent, parent_id, parent_source, pident, qcov, parent_coverage,
    evalue, bitscore), n_homolog_genes, n_self_hits.
    """
    out = pd.DataFrame(index=pd.Index(orfs["orf_id"], name="orf_id"))
    h, n_self = drop_self_hits(hits, orfs, spans)
    out["n_self_hits"] = n_self.reindex(out.index).fillna(0).astype(int)
    h = h[(h["evalue"] <= max_evalue) & (h["qcovhsp"] >= min_qcov)]
    h = h.sort_values(["qseqid", "bitscore", "pident"], ascending=[True, False, False])
    best = h.drop_duplicates("qseqid").set_index("qseqid")
    out["parent"] = best["parent"]
    out["parent_id"] = best["parent_id"]
    out["parent_source"] = best["source"]
    out["pident"] = best["pident"]
    out["qcov"] = best["qcovhsp"]
    out["parent_coverage"] = best["scovhsp"]
    out["evalue"] = best["evalue"]
    out["bitscore"] = best["bitscore"]
    out["n_homolog_genes"] = h.groupby("qseqid")["parent"].nunique().reindex(out.index).fillna(0).astype(int)
    return out


def run(orfs, gencode_exons, data_dir, workdir, diamond="diamond", threads=4):
    """End-to-end: reference set, DIAMOND search, homology table (index orf_id)."""
    d = os.path.join(data_dir, "homology")
    os.makedirs(workdir, exist_ok=True)
    ref_fa = os.path.join(workdir, "reference.faa")
    if not os.path.exists(ref_fa):
        write_fasta(reference_proteins(os.path.join(d, "gencode.v42.pc_translations.fa.gz"),
                                       [os.path.join(d, "UP000005640_9606.fasta.gz"),
                                        os.path.join(d, "UP000005640_9606_additional.fasta.gz")]), ref_fa)
    q = orfs[["orf_id", "aa_seq"]].dropna()
    query_fa = os.path.join(workdir, "orfs.faa")
    write_fasta(dict(zip(q["orf_id"], q["aa_seq"].str.rstrip("*"))), query_fa)
    hits_tsv = os.path.join(workdir, "hits.tsv")
    hits = read_hits(hits_tsv) if os.path.exists(hits_tsv) else run_diamond(
        query_fa, ref_fa, hits_tsv, workdir, diamond=diamond, threads=threads)
    return homology_table(orfs, hits, gene_spans(gencode_exons))
