"""
ENCODE HepG2 RNA-seq features of the mechanism classes (LOGIC.md, Step 4), from stranded signal bigWigs of
uniquely mapped reads (GRCh38), summed on the feature's strand over its exons minus the exons of other GENCODE
genes on the same strand (eclip.stranded_regions):

  nuclear enrichment  log2 (nucleus + eps) / (cytosol + eps), mean signal per base, replicates averaged
                      (Lecuyer CeFra-seq: nucleus ENCSR061SFU, cytosol ENCSR862HPO)

Features with too little signal in both conditions (mean signal per base below `min_signal`) are missing.
"""
import os

import numpy as np
import pandas as pd
import pyBigWig

from . import quant

ENCODE_DOWNLOAD = "https://www.encodeproject.org/files/{0}/@@download/{0}.bigWig"
# condition -> replicates -> (plus-strand, minus-strand) signal of unique reads
BIGWIGS = {
    "nucleus": [("ENCFF994NAN", "ENCFF935KNY"), ("ENCFF593IWR", "ENCFF614OZP")],
    "cytosol": [("ENCFF853VJV", "ENCFF704KFT"), ("ENCFF792BNC", "ENCFF194TSC")],
}
EPS = 0.01
MIN_SIGNAL = 0.05


def condition_signal(regions, strands, files, data_dir):
    """Mean signal per quantifiable base per feature, averaged over the replicates of one condition."""
    reg = regions.assign(strand=regions["feature_id"].map(strands)).dropna(subset=["strand"])
    n_bases = (reg["end1"] - reg["start0"]).groupby(reg["feature_id"]).sum()
    reps = []
    for plus, minus in files:
        total = []
        for acc, strand in ((plus, "+"), (minus, "-")):
            bw = pyBigWig.open(os.path.join(data_dir, f"{acc}.bigWig"))
            part = reg[reg["strand"] == strand].reset_index(drop=True)
            total.append(quant.sum_regions(bw, part).abs())
            bw.close()
        reps.append(pd.concat(total).groupby(level=0).sum() / n_bases)
    return pd.concat(reps, axis=1).mean(axis=1)


def log_ratio(a, b, eps=EPS, min_signal=MIN_SIGNAL):
    """log2 (a + eps) / (b + eps); missing where both a and b are below min_signal."""
    r = np.log2((a + eps) / (b + eps))
    return r.where((a >= min_signal) | (b >= min_signal))


def features(regions, strands, data_dir):
    """Per feature: mean signal per base in nucleus and cytosol (signal_) and nuclear_log2ratio."""
    sig = pd.DataFrame({c: condition_signal(regions, strands, f, data_dir) for c, f in BIGWIGS.items()})
    out = sig.add_prefix("signal_")
    out["nuclear_log2ratio"] = log_ratio(sig["nucleus"], sig["cytosol"])
    return out
