import numpy as np
import pandas as pd

from lincdisease import chromatin, eclip, eqtl


def test_bound_features_strand_and_overlap():
    regions = pd.DataFrame({"feature_id": ["a", "a", "b", "c"], "chrom": ["chr1"] * 4,
                            "start0": [100, 500, 100, 2000], "end1": [200, 600, 200, 2100]})
    strands = pd.Series({"a": "+", "b": "-", "c": "+"})
    peaks = pd.DataFrame({"chrom": ["chr1", "chr1", "chr1"], "start0": [150, 1000, 2100],
                          "end1": [180, 3000, 2150], "strand": ["+", "-", "+"]})
    # a: overlaps peak at 150; b: same position but other strand; c: long peak on '-' only, peak at 2100 is
    # adjacent (half-open) and does not overlap
    assert eclip.bound_features(regions, strands, peaks) == {"a"}
    peaks.loc[1, "strand"] = "+"
    assert eclip.bound_features(regions, strands, peaks) == {"a", "c"}


def test_localisation_groups():
    loc = pd.DataFrame(0, index=["N", "C", "B", "U"], columns=eclip.NUCLEAR + eclip.CYTOPLASMIC, dtype=float)
    loc.loc["N", "Speckles"] = 1
    loc.loc["C", "P bodies"] = 1
    loc.loc["B", ["Nuclei", "Cytoplasm"]] = 1
    loc.loc["U"] = np.nan
    assert eclip.localisation_group(loc).to_dict() == {"N": "nuclear", "C": "cytoplasmic", "B": "both", "U": ""}


def test_eclip_features_enrichment():
    regions = pd.DataFrame({"feature_id": ["x", "y"], "chrom": ["chr1", "chr1"], "start0": [0, 10_000],
                            "end1": [1000, 12_000]})
    strands = pd.Series({"x": "+", "y": "+"})
    pk = lambda s: pd.DataFrame({"chrom": ["chr1"], "start0": [s], "end1": [s + 50], "strand": ["+"]})
    sets = {"N1": pk(10), "N2": pk(20), "C1": pk(10_100), "O1": pk(10_200)}
    groups = pd.Series({"N1": "nuclear", "N2": "nuclear", "C1": "cytoplasmic", "O1": "both"})
    out = eclip.features(regions, strands, sets, groups)
    assert out.loc["x", "eclip_n_rbps"] == 2 and out.loc["y", "eclip_n_rbps"] == 2
    f = 0.5
    assert np.isclose(out.loc["x", "eclip_nuclear_enrichment"], np.log2(((2 + f) / 3) / f))
    assert out.loc["x", "eclip_nuclear_enrichment"] > 0 > out.loc["y", "eclip_nuclear_enrichment"]
    assert out.loc["y", "eclip_cytoplasmic_enrichment"] > 0
    assert np.isclose(out.loc["y", "eclip_rbps_per_kb"], 1.0)


def test_chromatin_scores():
    sig = pd.DataFrame({"H3K27ac": [10.0, 10.0, 0.5], "H3K4me1": [8.0, 1.0, 8.0], "H3K4me3": [1.0, 20.0, 1.0]},
                       index=["enh", "prom", "inactive"])
    s = chromatin.scores(sig)
    assert s.loc["enh", "enhancer_score"] > 0
    assert s.loc["enh", "enhancer_score"] > s.loc["prom", "enhancer_score"] > s.loc["inactive", "enhancer_score"]
    assert np.isclose(s.loc["enh", "enhancer_score"], np.log2(9 / 2))
    assert s.loc["inactive", "enhancer_score"] == chromatin.INACTIVE_SCORE


def test_eqtl_per_feature():
    egenes = pd.DataFrame({"gene_id": ["ENSG1.3", "ENSG2.1"], "gene_name": ["G1", "G2"],
                           "variant_id": ["chr1_1000_A_G_b38", "chr2_500_C_T_b38"],
                           "pval_nominal": [1e-10, 0.01], "slope": [0.5, 0.1], "qval": [1e-6, 0.4]})
    pairs = pd.DataFrame({"gene_id": ["ENSG1.3", "ENSG1.3"], "variant_id": ["chr1_1000_A_G_b38", "chr1_5000_G_C_b38"],
                          "pval_nominal": [1e-10, 1e-5], "slope": [0.5, 0.3]})
    gwas = pd.DataFrame({"chrom": ["chr1", "chr1"], "pos": [5000, 9000], "mlog10p": [12.0, 30.0]})
    loci = pd.DataFrame({"name_source_gene_ids": ["ENSG1.4", np.nan], "lncRNA_gene_ids": [np.nan, "ENSG9.1"]},
                        index=["L1", "L2"])
    links = eqtl.feature_genes(loci, ["G2"], egenes)
    out = eqtl.per_feature(links, egenes, pairs, gwas)
    assert set(out.index) == {"L1", "ctrl:G2"}            # L2's gene is not tested in GTEx
    assert out.loc["L1", "eqtl_egene"] and out.loc["L1", "eqtl_n_signif_pairs"] == 2
    assert out.loc["L1", "eqtl_gwas_shared"] and out.loc["L1", "eqtl_gwas_near"]
    assert out.loc["L1", "eqtl_gwas_mlog10p"] == 30.0      # 9000 is within 10 kb of the lead at 1000
    assert np.isclose(out.loc["L1", "eqtl_lead_mlog10p"], 10)
    assert not out.loc["ctrl:G2", "eqtl_egene"] and out.loc["ctrl:G2", "eqtl_gwas_mlog10p"] == 0


def test_stranded_regions_keep_antisense():
    lx = pd.DataFrame({"chrom": ["chr1"], "start0": [100], "end1": [1000], "strand": ["+"], "gene_id": ["L1"],
                       "gene_name": ["L1"], "gene_type": [""], "transcript_id": ["t"]})
    gc = pd.DataFrame({"chrom": ["chr1", "chr1"], "start0": [200, 600], "end1": [300, 700], "strand": ["-", "+"],
                       "gene_id": ["ENSGA.1", "ENSGS.1"], "gene_name": ["A", "S"], "gene_type": ["pc", "pc"],
                       "transcript_id": ["a", "s"]})
    loci = pd.DataFrame({"name_source_gene_ids": [np.nan], "lncRNA_gene_ids": [np.nan]}, index=["L1"])
    reg, strands = eclip.stranded_regions(lx, gc, loci, [])
    assert strands["L1"] == "+"
    assert (reg["end1"] - reg["start0"]).sum() == 800      # only the same-strand gene S is removed
