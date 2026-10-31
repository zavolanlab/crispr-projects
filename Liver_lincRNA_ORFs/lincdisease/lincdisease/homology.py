"""
Peptide homology and selection on the copy (LOGIC.md, Step 1).

Homology. Every ORF peptide is searched with DIAMOND blastp (--more-sensitive, default low-complexity masking)
against the GENCODE v42 protein-coding translations and the UniProt human reference proteome (canonical and
additional isoforms). Hits to the ORF's own locus (a GENCODE gene, or the gene named by a UniProt entry,
whose span overlaps the ORF) are discarded. A homolog is the best remaining hit (bitscore) with
E <= `max_evalue` and query coverage >= `min_qcov` %; reported as parent gene, protein id, identity (%),
query coverage (%), parent coverage (%) and the number of distinct homologous genes.

ORF integrity since duplication = parent coverage of the best hit (%): the fraction of the parent protein the
ORF still encodes before its stop codon / frameshift.

Selection on the copy. Ka/Ks (Nei-Gojobori 1986, Jukes-Cantor corrected) between the ORF coding sequence and the
CDS of the best GENCODE hit (whose CDS translates to the hit protein), on the codon alignment implied by the
DIAMOND protein alignment (columns with a gap, an ambiguous base or a stop codon are skipped). Ka/Ks is NaN
with fewer than `min_codons` aligned codons, without synonymous differences (dS = 0) or when p-distances are
saturated (>= 0.75).
"""
import gzip
import itertools
import os
import re
import subprocess

import numpy as np
import pandas as pd

BASES = "TCAG"
CODE = dict(zip((a + b + c for a in BASES for b in BASES for c in BASES),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

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


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def translate(nt):
    nt = nt.upper()
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def parse_blocks(blocks):
    """'s-e;s-e' (0-based, half-open) -> list of (start0, end1)."""
    return [tuple(int(x) for x in b.split("-")) for b in str(blocks).split(";") if b]


def orf_cds(orfs, genome_fa):
    """
    Coding sequence (stop codon included) of every ORF from its genomic blocks (`orf_blocks`, chrom, strand),
    reading the gzipped genome FASTA one chromosome at a time. Series orf_id -> sequence.
    """
    want = {c: g for c, g in orfs.groupby("chrom")}
    out = {}

    def extract(chrom, seq):
        for row in want[chrom].itertuples(index=False):
            s = "".join(seq[a:b] for a, b in parse_blocks(row.orf_blocks)).upper()
            out[row.orf_id] = revcomp(s) if row.strand == "-" else s

    opener = gzip.open if genome_fa.endswith(".gz") else open
    chrom, buf = None, []
    with opener(genome_fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if chrom in want:
                    extract(chrom, "".join(buf))
                chrom, buf = line[1:].split()[0], []
            elif chrom in want:
                buf.append(line.rstrip("\n"))
    if chrom in want:
        extract(chrom, "".join(buf))
    return pd.Series(out, dtype=object)


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


def gencode_cds(pc_transcripts):
    """Dict ENST id -> CDS (from the CDS:a-b field of the GENCODE pc_transcripts headers)."""
    out = {}
    for h, s in read_fasta(pc_transcripts).items():
        f = h.split("|")
        cds = [x for x in f if x.startswith("CDS:")]
        if cds:
            a, b = (int(x) for x in cds[0][4:].split("-"))
            out[f[0]] = s[a - 1:b]
    return out


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
    o = orfs.set_index("orf_id")[["chrom", "orf_start0", "orf_end1"]]
    m = hits.join(o, on="qseqid").merge(spans, left_on=["parent", "chrom"], right_on=["gene_name", "chrom"],
                                        how="left")
    is_self = (m["start0"] < m["orf_end1"]) & (m["end1"] > m["orf_start0"])
    self_rows = set(m.loc[is_self, "_row"])
    keep = ~hits["_row"].isin(self_rows)
    n_self = hits.loc[~keep].groupby("qseqid").size()
    return hits[keep], n_self


# ---------------------------------------------------------------- Ka/Ks

def _syn_sites(codon):
    s = 0.0
    aa = CODE[codon]
    for i in range(3):
        for b in BASES:
            if b != codon[i]:
                c = codon[:i] + b + codon[i + 1:]
                if CODE[c] == aa:
                    s += 1 / 3
    return s


SYN_SITES = {c: _syn_sites(c) for c in CODE if CODE[c] != "*"}


def _codon_diffs(c1, c2):
    """(synonymous, nonsynonymous) differences averaged over the mutational paths avoiding stop codons."""
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos:
        return 0.0, 0.0
    sd = nd = 0.0
    n_paths = 0
    for order in itertools.permutations(pos):
        cur, s, n, ok = c1, 0, 0, True
        for i in order:
            nxt = cur[:i] + c2[i] + cur[i + 1:]
            if CODE[nxt] == "*":
                ok = False
                break
            if CODE[nxt] == CODE[cur]:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            sd, nd, n_paths = sd + s, nd + n, n_paths + 1
    if n_paths == 0:
        return np.nan, np.nan
    return sd / n_paths, nd / n_paths


def codon_pairs(q_cds, s_cds, q_gapped, s_gapped, qstart, sstart):
    """Aligned codon pairs from a gapped protein alignment (1-based qstart/sstart of the alignment)."""
    i, j = qstart - 1, sstart - 1
    pairs = []
    for a, b in zip(q_gapped, s_gapped):
        if a != "-" and b != "-":
            c1, c2 = q_cds[3 * i:3 * i + 3].upper(), s_cds[3 * j:3 * j + 3].upper()
            if len(c1) == 3 and len(c2) == 3 and c1 in SYN_SITES and c2 in SYN_SITES:
                pairs.append((c1, c2))
        if a != "-":
            i += 1
        if b != "-":
            j += 1
    return pairs


def _jc(p):
    return -0.75 * np.log(1 - 4 * p / 3) if p < 0.75 else np.nan


def nei_gojobori(pairs, min_codons=20):
    """dict n_codons, S, N, pS, pN, dS, dN, kaks for aligned codon pairs."""
    res = dict(n_codons=len(pairs), S=np.nan, N=np.nan, pS=np.nan, pN=np.nan, dS=np.nan, dN=np.nan, kaks=np.nan)
    if len(pairs) < min_codons:
        return res
    S = sum((SYN_SITES[a] + SYN_SITES[b]) / 2 for a, b in pairs)
    N = 3 * len(pairs) - S
    d = [_codon_diffs(a, b) for a, b in pairs]
    sd = np.nansum([x[0] for x in d])
    nd = np.nansum([x[1] for x in d])
    pS, pN = sd / S, nd / N
    dS, dN = _jc(pS), _jc(pN)
    res.update(S=S, N=N, pS=pS, pN=pN, dS=dS, dN=dN)
    if np.isfinite(dS) and np.isfinite(dN) and dS > 0:
        res["kaks"] = dN / dS
    return res


# ---------------------------------------------------------------- per-ORF table

def homology_table(orfs, hits, spans, cds=None, parent_cds=None, max_evalue=1e-5, min_qcov=30.0,
                   min_codons=20):
    """
    One row per ORF: best non-self homolog (parent, parent_id, parent_source, pident, qcov, parent_coverage,
    evalue, bitscore), n_homolog_genes, n_self_hits, and Ka/Ks vs the best GENCODE hit with a usable CDS
    (kaks_parent_id, kaks_n_codons, pN, pS, dN, dS, kaks). `cds`: Series orf_id -> ORF CDS; `parent_cds`:
    dict ENST -> CDS.
    """
    out = pd.DataFrame(index=pd.Index(orfs["orf_id"], name="orf_id"))
    h = hits.reset_index(drop=True).assign(_row=lambda d: np.arange(len(d)))
    h, n_self = drop_self_hits(h, orfs, spans)
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

    cols = ["kaks_parent_id", "kaks_parent", "kaks_n_codons", "pN", "pS", "dN", "dS", "kaks"]
    for c in cols:
        out[c] = np.nan
    out["kaks_parent_id"] = out["kaks_parent_id"].astype(object)
    out["kaks_parent"] = out["kaks_parent"].astype(object)
    if cds is None or parent_cds is None:
        return out
    gc = h[h["source"] == "gencode"]
    for qid, g in gc.groupby("qseqid", sort=False):
        q = cds.get(qid)
        if not isinstance(q, str):
            continue
        for r in g.itertuples(index=False):
            s = parent_cds.get(r.parent_id)
            if s is None or not cds_matches(s, r.sseq_gapped, r.sstart):
                continue
            pairs = codon_pairs(q, s, r.qseq_gapped, r.sseq_gapped, r.qstart, r.sstart)
            ng = nei_gojobori(pairs, min_codons=min_codons)
            out.loc[qid, ["kaks_parent_id", "kaks_parent"]] = [r.parent_id, r.parent]
            out.loc[qid, ["kaks_n_codons", "pN", "pS", "dN", "dS", "kaks"]] = [
                ng["n_codons"], ng["pN"], ng["pS"], ng["dN"], ng["dS"], ng["kaks"]]
            break
    return out


def cds_matches(s_cds, s_gapped, sstart):
    """True if the translated parent CDS carries the aligned subject segment at `sstart` (1-based); residues
    masked by DIAMOND (X, low complexity) match anything."""
    seg = s_gapped.replace("-", "")
    prot = translate(s_cds)[sstart - 1:sstart - 1 + len(seg)]
    return len(prot) == len(seg) and all(a == b or b == "X" for a, b in zip(prot, seg))


def check_translation(orfs, cds):
    """Series orf_id -> True if the extracted CDS translates to the ORF peptide (aa_seq)."""
    aa = orfs.set_index("orf_id")["aa_seq"].astype(str).str.rstrip("*")
    tr = cds.reindex(aa.index).map(lambda s: translate(s).rstrip("*") if isinstance(s, str) else "")
    return tr.eq(aa)


def run(orfs, gencode_exons, data_dir, workdir, diamond="diamond", threads=4, genome_fa=None):
    """End-to-end: reference set, DIAMOND search, ORF CDS, homology + Ka/Ks table (index orf_id)."""
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
    cds = orf_cds(orfs, genome_fa or os.path.join(d, "GRCh38.primary_assembly.genome.fa.gz"))
    ok = check_translation(orfs, cds)
    res = homology_table(orfs, hits, gene_spans(gencode_exons), cds=cds[ok[ok].index],
                         parent_cds=gencode_cds(os.path.join(d, "gencode.v42.pc_transcripts.fa.gz")))
    res.insert(0, "cds_translates", ok.reindex(res.index))
    return res
