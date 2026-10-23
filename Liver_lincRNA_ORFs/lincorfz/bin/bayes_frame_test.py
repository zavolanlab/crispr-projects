#!/usr/bin/env python
"""
Bayesian test of in-frame vs out-of-frame Ribo-seq support of ORFs, with read-length-specific
priors on the frame distribution (frame_priors.py).

For an ORF and a dataset (one Ribo-seq sample, or the pooled samples), n_L = (n0, n1, n2) are the
numbers of covered P-site positions (default; reads stacking at one position are not independent
observations of the frame) or the P-site reads of aligned read length L in frames 0/1/2 of the ORF
(sense codons). The priors must be computed with the same convention. Models:

  H1        translation in the ORF frame:      n_L ~ DirMult(alpha_cds(L))
  H0 noise  untranslated background:           n_L ~ DirMult(alpha_null(L))      (3'UTR-like)
  H0 +1/+2  translation in another frame:      n_L ~ DirMult(alpha_cds(L) rotated by 1 or 2 frames)

The Dirichlet-multinomial marginalises over gene-to-gene variation of the frame distribution, so
each read length contributes evidence according to its own frame preservation: lengths whose
priors are close to uniform contribute ~0. Likelihoods are multiplied over read lengths (a single
scenario per ORF) and H0 is the mixture of its scenarios (weights --w-noise, and (1 - w_noise)/2
for each shifted frame). Reported per dataset:
  bayes_log10BF_<d>      log10 P(data | H1) / P(data | H0)
  bayes_posterior_<d>    posterior probability of H1 with prior --prior-h1
and for the pooled data additionally log10 BFs against each H0 scenario, the best alternative,
per-length counts / in-frame fractions / log10 BFs, the modal read length and the read lengths
that together carry >= --support-fraction of the ORF's P-sites.

Output: orf_frame_bayes.tsv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.special import gammaln, logsumexp

LN10 = np.log(10.0)


def dm_loglik(n, a):
    """Dirichlet-multinomial log-likelihood without the multinomial coefficient; n, a: (..., 3)."""
    A = a.sum(axis=-1)
    N = n.sum(axis=-1)
    return gammaln(A) - gammaln(N + A) + (gammaln(n + a) - gammaln(a)).sum(axis=-1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True,
                    help="orf_id, sample, length, f0, f1, f2 (reads), pos0, pos1, pos2 (covered positions)")
    ap.add_argument("--use", default="positions", choices=["positions", "reads"],
                    help="frame counts entering the test; must match the --position-cap of frame_priors.py")
    ap.add_argument("--priors", required=True)
    ap.add_argument("--prior-h1", type=float, default=0.5)
    ap.add_argument("--w-noise", type=float, default=0.5)
    ap.add_argument("--support-fraction", type=float, default=0.8)
    ap.add_argument("--min-log10bf-length", type=float, default=0.5)
    args = ap.parse_args()

    cnt = pd.read_csv(args.counts, sep="\t")
    pooled = cnt.groupby(["orf_id", "length"], as_index=False)[["f0", "f1", "f2", "pos0", "pos1", "pos2"]].sum() \
        .assign(sample="pooled")
    cnt = pd.concat([cnt, pooled], ignore_index=True)
    pri = pd.read_csv(args.priors, sep="\t").dropna(subset=["alpha_cds_0"])
    pri = pri.rename(columns={"dataset": "sample"})
    d = cnt.merge(pri[["sample", "length", "alpha_cds_0", "alpha_cds_1", "alpha_cds_2",
                       "alpha_null_0", "alpha_null_1", "alpha_null_2"]], on=["sample", "length"], how="inner")
    n = d[["pos0", "pos1", "pos2"] if args.use == "positions" else ["f0", "f1", "f2"]].to_numpy(float)
    a = d[["alpha_cds_0", "alpha_cds_1", "alpha_cds_2"]].to_numpy(float)
    a0 = d[["alpha_null_0", "alpha_null_1", "alpha_null_2"]].to_numpy(float)
    d["ll_h1"] = dm_loglik(n, a)
    d["ll_noise"] = dm_loglik(n, a0)
    d["ll_s1"] = dm_loglik(n, a[:, [2, 0, 1]])   # translation in frame +1 of the ORF
    d["ll_s2"] = dm_loglik(n, a[:, [1, 2, 0]])   # translation in frame +2 of the ORF
    lw = np.log([args.w_noise, (1 - args.w_noise) / 2, (1 - args.w_noise) / 2])

    def h0(ll_noise, ll_s1, ll_s2):
        return logsumexp(np.stack([ll_noise + lw[0], ll_s1 + lw[1], ll_s2 + lw[2]]), axis=0)

    # per read length (reported for the pooled data)
    d["log10bf_length"] = (d["ll_h1"] - h0(d["ll_noise"], d["ll_s1"], d["ll_s2"])) / LN10

    tot = d.groupby(["orf_id", "sample"])[["ll_h1", "ll_noise", "ll_s1", "ll_s2", "f0", "f1", "f2"]].sum()
    lh0 = h0(tot["ll_noise"].to_numpy(), tot["ll_s1"].to_numpy(), tot["ll_s2"].to_numpy())
    tot["log10BF"] = (tot["ll_h1"] - lh0) / LN10
    prior_odds = args.prior_h1 / (1 - args.prior_h1)
    tot["posterior"] = 1.0 / (1.0 + 1.0 / (prior_odds * np.power(10.0, tot["log10BF"].clip(-300, 300))))
    tot["log10BF_vs_noise"] = (tot["ll_h1"] - tot["ll_noise"]) / LN10
    tot["log10BF_vs_shifted_frame"] = (tot["ll_h1"] - logsumexp(
        np.stack([tot["ll_s1"], tot["ll_s2"]]), axis=0) + np.log(2)) / LN10
    alt = tot[["ll_noise", "ll_s1", "ll_s2"]].idxmax(axis=1).map(
        {"ll_noise": "untranslated", "ll_s1": "frame+1", "ll_s2": "frame+2"})
    tot["best_model"] = np.where(tot["ll_h1"] >= tot[["ll_noise", "ll_s1", "ll_s2"]].max(axis=1), "in_frame", alt)

    wide = tot[["log10BF", "posterior"]].unstack("sample")
    wide.columns = [f"bayes_{m}_{s}" for m, s in wide.columns]
    pt = tot.xs("pooled", level="sample")
    for c in ["log10BF_vs_noise", "log10BF_vs_shifted_frame", "best_model"]:
        wide[f"bayes_{c}_pooled"] = pt[c]

    # read-length composition of the pooled support
    pl = d.loc[d["sample"] == "pooled"].copy()
    pl["n"] = pl[["f0", "f1", "f2"]].sum(axis=1)
    pl = pl.loc[pl["n"] > 0].sort_values(["orf_id", "length"])
    rows = {}
    for oid, g in pl.groupby("orf_id"):
        n_tot = g["n"].sum()
        order = g.sort_values("n", ascending=False)
        cum = order["n"].cumsum() / n_tot
        k = int(np.searchsorted(cum.to_numpy(), args.support_fraction) + 1)
        sup = sorted(order["length"].iloc[:k].astype(int))
        good = g.loc[g["log10bf_length"] >= args.min_log10bf_length, "length"].astype(int).tolist()
        rows[oid] = {
            "modal_read_length": int(order["length"].iloc[0]),
            "support_read_lengths": ",".join(map(str, sup)),
            "support_read_length_range": f"{sup[0]}-{sup[-1]}",
            "read_lengths_with_inframe_evidence": ",".join(map(str, good)),
            "pooled_psites_by_length": ";".join(f"{int(L)}:{int(x)}" for L, x in zip(g["length"], g["n"])),
            "pooled_frac_inframe_by_length": ";".join(f"{int(L)}:{f0 / x:.2f}" for L, f0, x in
                                                      zip(g["length"], g["f0"], g["n"])),
            "pooled_log10BF_by_length": ";".join(f"{int(L)}:{b:.2f}" for L, b in zip(g["length"], g["log10bf_length"])),
        }
    out = wide.join(pd.DataFrame.from_dict(rows, orient="index"))
    out.index.name = "orf_id"
    out.reset_index().to_csv("orf_frame_bayes.tsv", sep="\t", index=False)
    print(out[["bayes_log10BF_pooled", "bayes_posterior_pooled"]].describe().to_string())
    print((out["bayes_log10BF_pooled"] >= 1).sum(), "ORFs with pooled log10 BF >= 1")


if __name__ == "__main__":
    main()
