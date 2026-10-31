"""
HepG2 chromatin at the dominant TSS of each feature (LOGIC.md, Step 4, M3 "enhancer-like chromatin").

Signal: ENCODE HepG2 histone ChIP-seq, fold change over control of the pooled replicates (GRCh38; BIGWIGS),
mean over TSS +- `flank` (1 kb). The dominant TSS is that of the isoform with the highest TPM in HepG2/Huh7
(`lincorfz` TPM_max) for lincRNA loci, and of the Ensembl canonical transcript for GENCODE controls.

  chrom_active    log2 ((H3K27ac + 1) / 2): > 0 when H3K27ac is enriched more than 1-fold above the control
  chrom_me1_me3   log2 ((H3K4me1 + 1) / (H3K4me3 + 1)): > 0 when H3K4me1 dominates H3K4me3
  enhancer_score  chrom_me1_me3 at active TSSs (H3K27ac fold change >= 2), INACTIVE_SCORE (-10, below any
                  observed log ratio) otherwise: highest for an active TSS with an enhancer-like (H3K4me1 >
                  H3K4me3) rather than a promoter-like signature, lowest for an inactive TSS
"""
import gzip
import os
import re

import numpy as np
import pandas as pd
import pyBigWig

ENCODE_DOWNLOAD = "https://www.encodeproject.org/files/{0}/@@download/{0}.bigWig"
# mark -> (experiment, pooled-replicate fold change over control)
BIGWIGS = {"H3K27ac": ("ENCSR000AMO", "ENCFF795ONN"), "H3K4me1": ("ENCSR000APV", "ENCFF576YVM"),
           "H3K4me3": ("ENCSR000AMP", "ENCFF284FVP")}
FLANK = 1000
MIN_H3K27AC = 2.0
INACTIVE_SCORE = -10.0


def dominant_tss(lincrna_tx, gencode_gtf, control_genes):
    """
    feature_id -> chrom/tss0/strand. `lincrna_tx`: lincRNA_transcripts.tsv (gene_id, start0, end1, strand,
    TPM_max); controls ('ctrl:<gene>') from the Ensembl canonical transcript of GENCODE.
    """
    t = lincrna_tx.sort_values("TPM_max", ascending=False).drop_duplicates("gene_id")
    loci = pd.DataFrame({"feature_id": t["gene_id"], "chrom": t["chrom"], "strand": t["strand"],
                         "tss0": np.where(t["strand"] == "+", t["start0"], t["end1"] - 1)})
    want = set(control_genes)
    rows = []
    opener = gzip.open if gencode_gtf.endswith(".gz") else open
    with opener(gencode_gtf, "rt") as fh:
        for ln in fh:
            if "\ttranscript\t" not in ln or "Ensembl_canonical" not in ln:
                continue
            f = ln.split("\t")
            g = re.search('gene_name "([^"]+)"', f[8]).group(1)
            if g in want:
                tss0 = int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1
                rows.append((f"ctrl:{g}", f[0], f[6], tss0))
    ctrl = pd.DataFrame(rows, columns=["feature_id", "chrom", "strand", "tss0"]).drop_duplicates("feature_id")
    return pd.concat([loci, ctrl], ignore_index=True).set_index("feature_id")


def mean_signal(bw, tss, flank=FLANK):
    chroms = bw.chroms()
    out = np.full(len(tss), np.nan)
    for i, (c, p) in enumerate(zip(tss["chrom"], tss["tss0"])):
        if c not in chroms:
            continue
        lo, hi = max(0, p - flank), min(chroms[c], p + flank + 1)
        v = bw.stats(c, lo, hi, type="mean", exact=True)[0]
        out[i] = 0.0 if v is None else v
    return pd.Series(out, index=tss.index)


def scores(sig):
    """`sig`: DataFrame with H3K27ac / H3K4me1 / H3K4me3 mean fold change."""
    out = sig.add_prefix("fc_")
    out["chrom_active"] = np.log2((sig["H3K27ac"] + 1) / 2)
    out["chrom_me1_me3"] = np.log2((sig["H3K4me1"] + 1) / (sig["H3K4me3"] + 1))
    active = sig["H3K27ac"] >= MIN_H3K27AC
    out["enhancer_score"] = out["chrom_me1_me3"].where(active, INACTIVE_SCORE).where(sig["H3K27ac"].notna())
    return out


def features(tss, data_dir, flank=FLANK):
    sig = {}
    for mark, (_, acc) in BIGWIGS.items():
        bw = pyBigWig.open(os.path.join(data_dir, f"{acc}.bigWig"))
        sig[mark] = mean_signal(bw, tss, flank)
        bw.close()
    out = scores(pd.DataFrame(sig))
    return tss.join(out)
