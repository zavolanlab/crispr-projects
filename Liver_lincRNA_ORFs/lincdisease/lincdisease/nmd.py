"""
NMD sensitivity of loci and controls (LOGIC.md, Step 4, M2): up-regulation upon inhibition of nonsense-mediated
decay, from recount3 coverage of the HeLa knock-downs of Colombo et al. 2017 (GSE86148, SRP083135), quantified like
TCGA/GTEx (same regions, `ncov`).

Samples: 6 scrambled-shRNA controls; NMD inhibited = UPF1 knock-down (3) and SMG6 + SMG7 double knock-down (3), the
two conditions in which NMD targets rise most. Single SMG6 / SMG7 knock-downs and the rescues are not used.

ncov is divided by median-of-ratios size factors of the recount3 gene sums (chronic.size_factors). Per feature,
nmd_log2fc = mean log2(scaled ncov + pseudocount) of the NMD-inhibited samples minus that of the controls
(> 0 = up when NMD is inhibited), as the mean of the two conditions' log2 fold changes; Welch t-test of the six
inhibited vs the six control samples, BH-FDR. Features detected (ncov > 0) in < `MIN_DETECTED` of the twelve
samples are missing: HeLa does not express most liver lincRNAs.
"""
import numpy as np
import pandas as pd
from scipy import stats

from . import quant
from .chronic import parse_attributes, size_factors
from .disease import bh

PROJECT = "SRP083135"
CONDITIONS = {"control": ("scramble", "none"), "UPF1_kd": ("UPF1", "none"), "SMG6_SMG7_kd": ("SMG6, SMG7", "none")}
MIN_DETECTED = 0.5


def manifest():
    """Samples of the control and NMD-inhibited conditions (quant.manifest format) with their condition."""
    man = quant.manifest("sra", PROJECT, ["sample_attributes", "sample_title"])
    a = parse_attributes(man["sample_attributes"]).set_index(man.index)
    kd, rescue = a["gene kd"].str.strip(), a["gene rescue"].str.strip()
    man["condition"] = pd.Series(np.nan, index=man.index, dtype=object)
    for name, (k, r) in CONDITIONS.items():
        man.loc[(kd == k) & (rescue == r), "condition"] = name
    return man.dropna(subset=["condition"]).reset_index(drop=True)


def sensitivity(ncov, samples, genes, pseudocount=0.05, min_detected=MIN_DETECTED):
    """
    ncov: features x samples (quant.normalize); samples: manifest rows with sample_id and condition; genes: AUC-scaled
    recount3 gene sums of protein-coding genes (genes x samples). Returns per feature nmd_log2fc, the two
    conditions' log2 fold changes, Welch p and BH-FDR.
    """
    cond = samples.set_index("sample_id")["condition"].reindex(ncov.columns)
    sf = size_factors(genes[ncov.columns])
    e = np.log2(ncov.div(sf, axis=1) + pseudocount)
    mean = {c: e.loc[:, cond == c].mean(axis=1) for c in CONDITIONS}
    out = pd.DataFrame(index=ncov.index)
    out["log2fc_UPF1_kd"] = mean["UPF1_kd"] - mean["control"]
    out["log2fc_SMG6_SMG7_kd"] = mean["SMG6_SMG7_kd"] - mean["control"]
    out["nmd_log2fc"] = out[["log2fc_UPF1_kd", "log2fc_SMG6_SMG7_kd"]].mean(axis=1)
    inhibited, control = e.loc[:, cond != "control"], e.loc[:, cond == "control"]
    out["nmd_p"] = stats.ttest_ind(inhibited, control, axis=1, equal_var=False).pvalue
    out["detected_fraction"] = (ncov > 0).mean(axis=1)
    out.loc[out["detected_fraction"] < min_detected, ["log2fc_UPF1_kd", "log2fc_SMG6_SMG7_kd", "nmd_log2fc",
                                                      "nmd_p"]] = np.nan
    out["nmd_fdr"] = bh(out["nmd_p"])
    return out
