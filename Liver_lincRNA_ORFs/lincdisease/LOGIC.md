# lincdisease: disease relevance and functional mechanisms of translated lincRNA ORFs

This file describes the current logic of the `lincdisease` module. It is kept as one clean, up-to-date
specification: when the logic changes, this file is edited in place (no change log; history is in git).

## Scope and inputs

* **Input**: a table of ORFs produced by `lincorfz` (same columns as `results/final/lincORFs.ranked.tsv`), taken as
  is. The current input is the manually curated `lincORFs.subselection.ranked.tsv.csv` of run `v2` (46 ORFs, 33 loci).
* **Columns of the input used here** (defined by `lincorfz`, `bin/integrate_and_rank.py` and the Sanity
  expression step): `orf_id`, `locus_id`, `locus_name`, `transcript_id`, ORF coordinates and `aa_seq`;
  `priority_score` and `rank` (the lincorfz ranking; candidates are ordered by it after the disease score),
  `ribo_evidence` (shown in the tables); `bayes_log10BF_pooled`, `frac_tx_psites_in_orf` and `frac_psites_unique`
  (translation, M1 and copy attribution); `MRL_pct_vs_mRNA` (polysome association, M3); `log2TPM_HepG2` and
  `log2TPM_Huh7` (cell-line expression, Sanity estimates without pseudocount, so negative values occur; D6, M4); and
  the lincorfz QC flags listed in Step 1.
* **Background**: all expressed lincRNA loci and their ORFs from the same `lincorfz` run. Every measure is also
  computed for the background, and scores are expressed as percentiles against it (0.5 = a typical expressed lincRNA).
* **Disease focus**: liver disease, with hepatocellular carcinoma (HCC) as the main but not the only focus
  (MASLD/MASH and fibrosis, viral hepatitis, liver-function traits).
* **Two levels**: disease evidence attaches mostly to the lincRNA **locus**; mechanism evidence separates what the
  **RNA** does from what the **peptide** does. Locus and ORF scores are reported side by side and never multiplied.

## General rules

1. Each component is a statistic converted to an empirical percentile among the background. Signed statistics
   are ranked on magnitude, so a locus strongly down in disease scores as high as one strongly up; the sign is kept
   separately as a direction (`*_direction`). Where the background is ORFs (mechanism,
   actionability), each ORF is weighted by 1 / number of ORFs of its locus, so that every locus counts once:
   locus-level features are copied to all ORFs of a locus, and loci with many ORFs would otherwise dominate.
2. Missing data never counts as zero. A composite score is S = sum_k w_k q_k / sum_k w_k, with q_k the percentile
   of component k, or 0.5 if it is missing; the weights w_k are fixed by hand and listed with each score (disease
   score below, D1 and D6 in Step 3, mechanism classes in Step 4). Its `coverage` is C = sum_k w_k c_k / sum_k w_k,
   with c_k = 1 if component k is observed, 0 if missing, and its own coverage if it is itself a composite (D1,
   D6: a locus with D1 HCC but no D1 chronic has c_D1 = 0.7). Then S = C * S_obs + (1 - C) * 0.5, where S_obs is
   the weighted mean of the observed percentiles: the observed evidence shrunk toward 0.5 in proportion to the
   missing weight.
3. Quality checks produce flags and a confidence tier; they do not remove ORFs and do not change the scores.
4. Default weights are fixed in advance. Fitted weights (see "Weight fitting") are reported as an alternative,
   never silently substituted.

## Overview of the scores

| ID | Name | Statistic (per locus) | Signed | Centred | Weight in `disease_score` | Use |
|:---|:---|:---|:---|:---|:---|:---|
| D1 HCC | Tumour vs adjacent normal | paired rank-biserial r, TCGA-LIHC | yes | yes | 0.175 (0.25 x 0.7) | scored |
| D1 chronic | Chronic liver disease | Stouffer z over three biopsy cohorts | yes | yes | 0.075 (0.25 x 0.3) | scored |
| D1 replication | Four independent HCC cohorts | combined paired r, `hcc_replicated` | yes | compared to the median | - | reported |
| D2 | Prognosis | mean Wald z of OS and PFI | yes | yes | 0.25 | scored |
| D3 | Somatic copy number | dosage response beyond the neighbouring genes x (gain - loss frequency) | yes | no | 0.20 | scored |
| D4 | Germline genetics | max -log10 p of liver-trait GWAS hits within 50 kb | no | no | 0.10 | scored (eQTLs reported) |
| D5 | Prior knowledge | log2(1 + liver papers) | no | no | 0.10 | scored (`disease_score_noD5` without it) |
| D6 | Liver context | composite of liver expression, liver specificity, HepG2/Huh7 expression | no | no | 0.10 | scored |
| M1-M6 | Mechanism classes (per ORF) | mean feature percentile and rank-based z | no | - | - | separate from the disease score |

Signed: ranked on magnitude, sign kept as the direction. Centred: the background median is subtracted before
ranking (Step 3). `hcc_score` and `chronic_score` reweight the same components (Step 3).

## Glossary

* **D1-D6**: disease components (Step 3); **M1-M6**: mechanism classes (Step 4): M1 encoded peptide, M2
  translation-coupled RNA regulation, M3 cis-regulation, M4 nuclear trans-acting RNA, M5 cytoplasmic trans-acting
  RNA, M6 small-RNA host.
* **ncov**: normalized coverage of a locus in a sample (Step 2), an RPKM-like unit.
* **Total coverage** (recount3 `bc_auc`): the summed per-base coverage of a sample, used to scale ncov. Not to be
  confused with **AUC** in the statistical sense, the area under the ROC curve (`auc_vs_gtex`; the cross-validated
  AUC of the weight fit).
* **TPM / log2TPM**: transcripts per million; the lincorfz cell-line values are Sanity log2 estimates without
  pseudocount.
* **P-sites**: the ribosomal P-site positions inferred from Ribo-seq reads; **MRL**: mean ribosome load (polysome
  profiling, lincorfz).
* **NMD**: nonsense-mediated mRNA decay. **TSS**: transcription start site. **pgRNA**: processed pseudogene RNA.
* **CeFra-seq**: RNA-seq of subcellular fractions (ENCODE HepG2 nucleus / cytosol). **IDR**: irreproducible
  discovery rate, ENCODE's replicate-consistency filter for eCLIP peaks.
* **NAS**: NAFLD activity score; **OS / PFI**: overall survival / progression-free interval; **TCGA-CDR**: the TCGA
  Pan-Cancer Clinical Data Resource (curated endpoints).
* **eGene**: a gene with a significant eQTL; **LD**: linkage disequilibrium.
* **BH-FDR**: Benjamini-Hochberg false discovery rate.
* **Rank-biserial correlation**: the effect size of a rank test, (favourable - unfavourable rank sum) / total, in
  [-1, 1]. **Stouffer z**: the weighted sum of per-cohort z-scores divided by the root of the summed squared
  weights.
* **Boldanova 2017** is the HCV cohort; its columns are named `hcv_*`.

## Step 1. Quality-check layer (per ORF)

Flags:

* overlap of the ORF's coding exons (not its introns) with segmental duplications (UCSC `genomicSuperDups`, with
  identity `fracMatch`) and with GRCh38 regions that are falsely duplicated relative to T2T-CHM13;
* peptide homology (DIAMOND blastp `--more-sensitive` against GENCODE v42 protein-coding translations and the
  UniProt human reference proteome with isoforms; hits to genes overlapping the ORF are discarded; homolog = best hit
  with E <= 1e-5 and query coverage >= 30%): parent gene, identity, query coverage, parent coverage (ORF integrity
  since duplication; reported, not scored, as a genuine short ORF of a pseudogene covers little of its parent),
  number of homologous genes;
* flags already computed by `lincorfz`: antisense overlap with coding exons, small-RNA overlap,
  `resolved_multimapper_dominated`, possible 5'/3' extension of a coding gene;
* `tcga_quantifiable`: TCGA and GTEx RNA-seq is unstranded, so for loci antisense to expressed genes only exonic
  positions not overlapping the sense gene are quantified; loci with too few such positions are flagged.

**Copy-attribution confidence** (0-1), for ORFs in duplicated regions or with homology to annotated proteins:

1. copy-specific translation: fraction of P-sites from uniquely mapped reads;
2. copy-specific peptides: identified peptides (proteomics, immunopeptidomics) covering residues that differ from
   the paralogs (1 if present; missing otherwise, as absence of peptides is not evidence against the copy).
   Implemented as a peptide absent from canonical proteins and contaminants that maps to this lincORF only;
3. alignment attribution: age of the duplication, (0.98 - identity) / 0.08 clipped to 0-1, i.e. 0 at >= 98% and 1
   at <= 90% identity (identity = segmental
   duplication `fracMatch`, else protein identity to the homolog). For old duplications the multi-species
   alignment likely contains true orthologs of the copy; for very recent ones (> 98% identity) other species
   probably carry the parent.

The confidence is the mean of the available parts (1 for ORFs outside duplications and without homologs).
Parts 1-2 are the evidence that the copy itself is translated; part 3 only says how far the conservation
alignment belongs to the copy.

Conservation is **not** discarded in duplicated regions: if the copy is expressed and translated, conservation
indicates evolutionary support for its function. The M1 features that a paralogous copy can supply (translation
strength, P-site share, sequence features, PhyloCSF, PhyloP contrast) enter weighted by the copy-attribution
confidence (Step 4).

**Confidence tier**: `clean`; `caution` (antisense overlap with coding exons, small-RNA overlap,
multimapper-dominated, possible coding-gene extension, not quantifiable in TCGA/GTEx); `expressed_paralog`
(duplicated or homologous, with copy-specific evidence); `likely_artefact` (duplicated or homologous, and the
translation evidence cannot be attributed to the copy: < 50% unique P-sites or multimapper-dominated, and no
copy-specific peptide).

## Step 2. Quantification in patient and tissue data

Many loci are de novo assemblies without a GENCODE gene, so gene-level public tables cannot be used. Coverage is
read from recount3 base-level coverage bigWigs, per locus (disease evidence attaches to the locus), over its
**quantifiable positions**: the union of the exons of its isoforms minus the exons of every GENCODE gene on
either strand that is not the locus itself (public RNA-seq is unstranded). Controls are measured the same way on
the union of their GENCODE exons. Loci with < 100 quantifiable positions are flagged `tcga_quantifiable = False`.
Samples:

* TCGA-LIHC (all tumours and adjacent normals), and four independent HCC cohorts of tumours and adjacent normals in
  recount3 (`quant_hcc/`, D1 HCC replication),
* GTEx (all normal liver samples; 10 random samples of every other tissue),
* recount3 SRA liver disease cohorts (`quant_sra/`): Hoang 2019 (GSE130970, 78 NAFLD biopsies with fibrosis stage and
  NAFLD activity score), Suppli 2019 (GSE126848, 57: normal-weight, obese, NAFL, NASH), Boldanova 2017 (GSE84346,
  untreated biopsies, HCV vs non-HCV); phenotypes from the recount3 SRA `sample_attributes`. Govaere 2020
  (GSE135251, 206 NAFLD biopsies with fibrosis stage and NAS, 10 controls; not in recount3): STAR alignment
  (GENCODE v42), loci counted in the CRAMs over the same regions (samtools bedcov), gene counts and technical
  covariates (input reads, % uniquely mapped, exonic fraction = reads assigned to genes / uniquely mapped reads) from
  STAR, phenotypes from the GEO sample characteristics.

Values (`ncov`): summed coverage / number of quantifiable positions / the sample's total base coverage
(recount3 `bc_auc.all_reads_all_bases`) x 1e9, divided by a median-of-ratios **size factor** of the sample
(protein-coding recount3 gene sums detected in >= 90% of the samples, the most abundant 1% left out; geometric mean
1). Total coverage alone would let the few hepatocyte transcripts (ALB, APOA1, ...) that dominate liver RNA set the
scale; their share falls in tumours and injured liver, which would lift every other transcript. Size factors are
computed over the samples that are compared with each other: TCGA-LIHC and GTEx together (`scaling` step,
`quant/size_factors.tsv`: tumours, adjacent normals and the GTEx tissues are compared in D1, D6 and tumour
specificity), and per cohort for the HCC replication, chronic liver disease and NMD cohorts (Govaere: summed
coverage and STAR gene counts; its statistics are within-cohort ranks). Statistics use log2(ncov + pseudocount). The recount3 gene sums and the GTEx eQTLs use GENCODE v26 ids, the
regions GENCODE v42; genes are matched on unversioned Ensembl ids.
Loci flagged `tcga_quantifiable = False` are not tested in any expression-based statistic and stay out of its
multiple-testing (BH) family. Every sums directory keeps the regions it was computed on (`regions.tsv`); reading
sums whose regions differ from the current `quant/regions.tsv`, or that lack a feature, is an error.

## Step 3. Disease components (per locus)

| | Component | Measure |
|:---|:---|:---|
| D1 | Disease-associated expression | HCC (D1 HCC): paired tumour vs adjacent normal, TCGA-LIHC. Chronic liver disease (D1 chronic): association with fibrosis or infection in three biopsy cohorts. HCC carries 0.7 and chronic liver disease 0.3 of D1. Details below. |
| D2 | Prognosis | Cox models of overall survival and progression-free interval in TCGA-LIHC (TCGA-CDR endpoints; one primary tumour per patient) on standardized expression, adjusted for age, stage (I-II vs III-IV) and the tumour's log size factor (less differentiated tumours lose hepatocyte transcripts, and this composition shift itself predicts survival); loci expressed in < 20% of tumours are not tested; statistic = mean Wald z of OS and PFI (HR per SD and BH-FDR reported; OS and PFI are correlated endpoints, so D2 is close to one endpoint) |
| D3 | Somatic copy number | TCGA-LIHC primary tumours (GDC GRCh38 masked copy-number segments from the UCSC Xena GDC hub): segment value (log2 ratio) at the midpoint of the locus span; Spearman rho of CN vs expression, partial on tumour purity (ABSOLUTE, TCGA PanCanAtlas; rank residuals after linear regression on purity rank, tumours with a purity call), as non-tumour cells shift expression and CN amplitude together (unadjusted rho reported); gain / loss frequency (> +0.3 / < -0.3). HCC copy-number changes span whole arms (1q, 8q gains), where most expressed genes follow the copy number to some degree; the dosage response is therefore measured relative to the neighbourhood: the same purity-adjusted rho for the protein-coding genes (recount3 gene sums of the same tumours, scaled by the size factors, detected in >= 20% of them), `cn_rho_neighbours` = median over the 20 nearest on the same chromosome, `cn_rho_excess` = rho_CN - `cn_rho_neighbours`. Statistic = max(`cn_rho_excess`, 0) x (gain - loss frequency), and 0 where the CN-expression correlation is not significant (`cn_fdr` >= 0.05): > 0 when expression follows a frequent gain more closely than the genes around it, < 0 when it follows a frequent loss. |
| D4 | Germline genetics | GWAS Catalog associations with whitelisted liver traits (EFO labels of HCC, MASLD/MASH, liver fibrosis, cirrhosis, alcoholic liver disease, viral hepatitis, ...; ALT, AST, GGT) within the locus span +- 50 kb (1-based positions start + 1 - 50 kb to end + 50 kb; an LD-window proxy); statistic = largest -log10 p (0 without association), reported per trait class (disease / enzyme). GTEx v8 liver eQTLs for loci with GENCODE genes (v26 ids without version): eGene q < 0.05, lead variant, and a link to a whitelisted GWAS hit if a significant eQTL variant is the hit or the lead eQTL lies within 10 kb (distance proxy, no LD panel); reported (`D4.eqtl.tsv`, `D4_eqtl_supported`), not scored, as only annotated loci can have GTEx eQTLs |
| D5 | Prior knowledge | Europe PMC, title or abstract: exact symbol or an HGNC previous symbol / alias (up to 5), AND an RNA term (guards against symbols that are words), AND a liver term; counts of all, liver, review and retracted papers; statistic = log2(1 + liver papers - retracted). Usable aliases: >= 5 characters with a digit, not another gene's approved symbol, not listed by HGNC for more than one gene, not the stem of another gene's symbol (EEF1B, the protein name of EEF1B2), not a generic token (a Greek letter or 1-3 letters with a number, e.g. BETA-3); words like ALIEN or short tokens like LEU2 would collect unrelated papers. Europe PMC splits phrases at hyphens, so a `<gene>-AS/-IT/-DT[n]` name is dropped from the query when it has >= 0.5 times the liver papers of `<gene>` (which has >= 20); a locus with no usable name left has missing counts; dropped names are reported (`names_dropped`). Exact-phrase matching keeps papers on LINC01089 from being credited to LINC01090. Loci named only by an Ensembl, StringTie or NONCODE (NONHSAG) id count 0. |
| D6 | Liver context | Equal-weight composite of the percentiles of: median expression in GTEx liver; liver specificity (liver median minus the largest median of the other GTEx tissues); max(log2 TPM HepG2, Huh7) from `lincorfz` |

### D1 HCC

Matched-pairs rank-biserial correlation r of TCGA-LIHC primary tumours vs adjacent normals of the same patients
(Wilcoxon signed-rank test; > 0 = up in tumour). Zero differences, common for lowly expressed loci, are kept in the
ranking (Pratt): all |differences| are ranked, zeros included, and r = (positive - negative rank sum) / sum of all
ranks, so a locus undetected in most patients cannot reach an extreme r from its few detected pairs. Loci detected
(ncov > 0) in < 20% of the paired samples are not tested. log2FC and BH-FDR are reported. Tumours vs GTEx liver
(Mann-Whitney; AUC = (r + 1) / 2) is reported, not scored: TCGA adjacent normals and GTEx liver, both non-tumour
liver, already differ (|r| > 0.3) for most loci (study batch), and the paired design is the only one free of it.

**Replication** (reported, not scored): the same matched-pairs statistic in four independent recount3 cohorts of
tumours and adjacent non-tumour liver of the same patients: GSE77314 (49 pairs), GSE77509 (20 HBV-HCC with portal
vein tumour thrombus; runs of a sample averaged, thrombi not used), GSE94660 (21 HBV-HCC), GSE56545 (14). A feature
detected in < 20% of a cohort's samples is not tested there. Combined r = pair-weighted mean, Stouffer p (weights
sqrt(pairs)), BH-FDR. `hcc_replicated` = combined FDR < 0.05 and the combined r and D1 HCC on the same side of
their background medians (missing if no cohort tested the locus).

### D1 chronic liver disease

Per cohort, a partial rank correlation with the phenotype (> 0 = higher with disease): ranks of expression and
phenotype residualized on the ranks of log input reads, % uniquely mapped reads and exonic fraction, t test with
n - 2 - k df. Phenotypes: Hoang and Govaere fibrosis stage (NAFLD activity score reported), Boldanova HCV vs
non-HCV, Suppli NAFL + NASH vs normal-weight + obese. A feature detected (ncov > 0) in < 20% of a cohort's samples
is not tested there. Statistic `D1_chronic` = Stouffer z of the scored cohorts (Hoang, Govaere, Boldanova): the z of
each cohort's two-sided p with the sign of its correlation, weights sqrt(n tested), so that each cohort counts by
its size; BH-FDR of its p. `chronic_r` (reported) = the n-weighted mean of the cohort correlations.

Suppli is reported, not scored: its cases differ technically from its controls (fewer uniquely mapped reads, lower
exonic fraction; the 15 samples with the lowest exonic fraction are all cases), and after rescaling and adjustment
40% of loci keep p < 0.01 (median r +0.14; Hoang 9%, 0.00).

Parent genes: a locus whose expression follows a parent gene (homology step; median Spearman over cohorts >= 0.5)
that is itself associated with chronic disease (FDR < 0.05, same sign) is `parent_driven`, and its cohort
statistics, scored and cell-intrinsic, are recomputed with the parent's expression as an extra covariate
(processed pseudogenes pick up reads of parents such as KRT8, KRT18, YWHAZ that rise in injured liver); the parent
genes' own association is BH-corrected over the parent genes tested.

Reported, not scored: `D1_chronic_cell`, the same statistic additionally adjusted for liver cell-type marker
scores (hepatocyte, stellate, cholangiocyte, macrophage, lymphocyte, endothelial; MacParland et al. 2018 markers,
COL1A1, KRT8, KRT18 left out), i.e. the part not explained by a change in cell mix.

### Centring

Across the background most loci are somewhat higher in tumours than in their adjacent normals, somewhat lower in
chronic disease and somewhat higher with worse survival. A locus that only follows such a shift is not
specifically disease-associated, so D1 HCC, D1 chronic and D2 are ranked, and give their direction, after
subtracting their median over the background loci. The `stat_` columns keep the uncentred statistics.

### Disease score

**Disease score** = weighted mean of component percentiles, shrunk toward 0.5 for missing components.
Default weights: D1 0.25, D2 0.25, D3 0.20, D4 0.10, D5 0.10, D6 0.10. Reported with:

* `disease_direction`: from the signs of the informative signed components. A component (D1 HCC, D1 chronic,
  D2 or D3, each centred as above where it applies) is informative if its own percentile is >= 0.75;
  oncogene-like if all informative signs are positive, suppressor-like if all are negative, inconsistent if both,
  none if no component is informative;
* `disease_coverage`, the confidence tier, and `disease_score_noD5` (without prior knowledge; the other
  weights rescaled to sum 1). Candidates are ordered by `disease_score_noD5`, so that the ranking does not favour loci that
  are already well studied.

**HCC score and chronic liver disease score**, reported next to the disease score (same percentiles, weights
fixed in advance; the disease score is mostly an HCC score, with chronic disease expression at 7.5%):

| Score | Components (weight) |
|---|---|
| `hcc_score` | D1 HCC (0.25), D2 (0.25), D3 (0.20), D4 HCC GWAS (0.10), D5 liver cancer papers (0.10), D6 (0.10) |
| `chronic_score` | D1 chronic (0.50), D4 chronic liver disease GWAS and liver enzymes (0.20), D5 chronic liver disease papers (0.15), D6 (0.15) |

D4 split: HCC = the "hepatocellular carcinoma" trait; chronic = the other whitelisted liver disease traits, and the
liver enzymes (max -log10 p). D5 split: the liver papers that also match a cancer term (cancer, carcinoma, tumour,
HCC, hepatoma, hepatocellular, oncogen*, metasta*) or a chronic disease term (NAFLD, MASLD, NASH, MASH, steatosis,
steatohepatitis, fatty liver, fibrosis, cirrhosis, hepatitis, liver injury / failure, cholestasis); a paper can
count for both; statistic log2(1 + papers). Directions: `hcc_direction` from D1 HCC, D2 and D3 (as
`disease_direction`); `chronic_direction` from D1 chronic (up / down with disease, percentile >= 0.75).

**Actionability score** (per ORF, separate from the disease score): equal-weight composite of the percentiles
(among all ORFs, locus-weighted) of tumour-specific expression of the locus (median in TCGA-LIHC primary tumours minus the
largest median over GTEx tissues except testis), predicted HLA class I presentation (MHCflurry 2 presentation
percentile <= 2, distinct 8-11-mer binders over HepG2 A*02:01, A*24:02, B*35:14, B*51:01, C*04:01, C*16:02, Huh7
A*11:01, B*54:01, C*01:02 (Cellosaurus) and 16 common alleles, each scored on its own; peptides that also occur in
an annotated protein, GENCODE v42 translations or UniProt human with isoforms, are self-peptides and are left out),
and the number of lincORF-specific HLA-I immunopeptides (`lincorfz` notebook section 7.3).

## Step 4. Mechanism classes (per ORF, with locus-level features)

Features are ranked against all ORFs of the run, each ORF weighted by 1 / number of ORFs of its locus; locus-level
features are copied to every ORF of the locus.
ENCODE HepG2 features use stranded data, so they are measured on the exons of the feature minus the exons of other
GENCODE genes on the same strand only (`quant/regions.tsv` also removes antisense exons, which unstranded RNA-seq
needs but stranded data do not). Fractionation uses stranded signal bigWigs of uniquely mapped reads (two
replicates, averaged), as mean signal per position on the locus strand: nuclear log2 ratio = log2 (nucleus +
0.01) / (cytosol + 0.01) from CeFra-seq; missing if both conditions have < 0.05 signal per base. Its negative,
cytoplasmic enrichment, is a feature of M5 only: it defines a cytoplasmic RNA, whereas cis-acting (M3) and nuclear
trans-acting (M4) RNAs are both nuclear, so it cannot tell those two apart.

Ranking adjustments:

* Copy attribution: the M1 features that a paralogous copy can supply as well as this one (translation strength,
  P-site share, sequence features, PhyloCSF, PhyloP contrast) enter as 0.5 + (percentile - 0.5) x copy-attribution
  confidence (`pct_<feature>_weighted`).
* ORF length: features that rise with ORF length (translation strength, P-site share, sequence features, PhyloCSF;
  locus-weighted Spearman 0.19-0.56 with length) are ranked among the ORFs of the same length bin (20 quantile bins
  of the background ORF lengths; `length_bin`), so that M1 does not reward length as such; length itself is not a
  feature.
* HepG2 expression: features that rise or fall with HepG2 expression (PhyloCSF, NMD sensitivity, low polysome
  association; |Spearman| 0.34-0.44 with log2 TPM HepG2) are ranked among the ORFs of the same expression bin (5
  quantile bins; `expression_bin`; ORFs without HepG2 expression form a bin of their own). PhyloCSF, in both lists,
  is ranked in joint length x expression bins. M4 abundance is itself the expression and is ranked as it is.
HepG2 eCLIP (M4, M5): reproducible (IDR) peaks of 105 RBPs (ENCODE), log2 fold enrichment >= 3 and -log10 p >= 3,
on the locus strand. RBP groups by HepG2 immunofluorescence (Van Nostrand et al. 2020, Suppl. Data 1): nuclear =
only nuclear compartments (17 RBPs), cytoplasmic = only cytoplasmic compartments (15); RBPs seen in both (61) or not
imaged (12) count only in the total. Statistic per group: log2[((n_group + f) / (n_bound + 1)) / f], f = the group's
share of the 105 RBPs (over-representation among the bound RBPs, shrunk toward 0 by one pseudo-RBP: 0 without
binding, further from 0 the more RBPs support it); loci expressed in HepG2 (TPM >= 1) without a bound RBP have 0, no
evidence either way; loci not expressed in HepG2 are missing.
Each class has an absolute evidence score (mean of its feature percentiles, missing features counting 0.5) and a
rank-based standardized score z, per class: the standard normal quantile of the percentile of the class score among
the ORFs of the run that have data for the class (locus-weighted). Classes are compared on z, because their raw
scores have different spreads (few or binary features give wide ones); the rank transform maps the top of every
class to the same quantiles, so that a class with one rare binary feature cannot win by its spread alone. An ORF
without data for a class (coverage 0) has no z for it, and that class does not compete. Classes are not mutually
exclusive; a class is called dominant only if its z is > 0 (above the median) and exceeds the next class's z by
>= 0.5.

| Class | Features |
|:---|:---|
| M1 Encoded peptide | Translation strength (pooled log10 Bayes factor) and fraction of the transcript's P-sites in the ORF (`lincorfz`); identified peptides (proteomics peptides absent from canonical proteins, plus lincORF-specific immunopeptides); independent translation evidence: same stop codon as an ORF of the GENCODE Ribo-seq ORF catalog (phase 2 comprehensive, v45), +1 if it matches a phase 1 ORF with peptide support; PhyloCSF; ORF vs non-ORF exon conservation: mean PhyloP100way over the ORF minus the mean over the exonic positions of its transcript outside the ORF, both restricted to the locus' quantifiable positions (Step 2: exons of other GENCODE genes on either strand removed, as PhyloP is strand-blind and antisense coding exons or miRNA hairpins would lend their conservation to the ORF; missing with < 30 ORF or < 30 other positions); sequence features = 0.5 x [Kyte-Doolittle TM segment (19-aa window mean >= 1.6), signal-peptide-like N-terminus (heuristic: hydrophobic 8-mer, mean >= 2.0 without D/E/K/R, starting at residues 2-20, small residues at -3/-1 of a cleavage site <= 40; open predictors only, SignalP/TMHMM are licensed), or Pfam-A domain (hmmsearch, gathering thresholds)] + 0.5 x (1 - fraction of residues disordered, metapredict >= 0.5). Published CRISPR knock-out screens, matched by the same stop codon: (a) Prensner et al. 2021 (553 non-canonical ORFs; primary screen in 8 cell lines including HepG2, secondary tiling screen of the hits): 2 = viability phenotype specific to the ORF (secondary screen: specific to the ORF, or to the ORF and one other transcript region), 1 = phenotype not resolved by the secondary screen, 0 = no phenotype or a phenotype attributed to the transcript, the genomic region or an adjacent coding gene; the HepG2 median sgRNA log2 fold change is reported. (b) Chen et al. 2020 (~2,300 ORFs, K562 and iPSC; coordinates lifted from hg19 with the UCSC chain): 1 = hit in either cell line, 0 = screened without a hit. Their hit calls are not in the supplementary tables, and their stated null (non-targeting sgRNAs) makes nearly every ORF significant, because Cas9 cutting alone costs a phenotype of about -0.015 (their estimate); hits are re-called as BH-FDR < 0.05 in a one-sided Mann-Whitney test of the ORF's sgRNA phenotypes against the non-targeting sgRNAs shifted by -0.015, and a top-3 sgRNA phenotype <= -0.02 (808 hits in either cell line, 146 in both; the paper reports 401 in both, so the re-call is stricter). Feature: Prensner's secondary-screen verdict where it exists (the only test of ORF specificity), otherwise the larger of the two values; missing if neither screen tested the ORF. Translation strength, P-site share, sequence features, PhyloCSF and PhyloP contrast are weighted by copy attribution (see above). |
| M2 Translation-coupled RNA regulation | Stop codon > 50 nt upstream of the last exon junction (NMD substrate); NMD sensitivity: up-regulation upon NMD inhibition in HeLa (Colombo et al. 2017, GSE86148 / recount3 SRP083135; quantified like TCGA/GTEx over the locus' quantifiable positions, ncov scaled by median-of-ratios size factors of the protein-coding gene sums; mean of the log2 fold changes of the UPF1 knock-down and the SMG6 + SMG7 double knock-down, 3 replicates each, vs 6 scrambled controls; t test of this contrast in a linear model with the three conditions (residual variance pooled within conditions, 9 df), BH-FDR; missing if the locus is detected in < 50% of the 12 samples). The feature is the fold change shrunk toward 0 by its standard error (empirical Bayes, normal prior): raw x tau^2 / (tau^2 + se^2), tau^2 = variance of the raw fold changes minus their median se^2 (`nmd_log2fc`; `nmd_log2fc_raw` and `nmd_se` reported), so that noisy, lowly covered loci do not reach the top by chance |
| M3 Cis-regulation of a neighbour (RNA or DNA element) | Divergent (head-to-head, TSSs within 5 kb) or antisense (exon overlap) position next to a protein-coding gene; co-expression with the nearest protein-coding gene (TSS distance <= 100 kb, either strand; recount3 gene sums) across TCGA-LIHC and GTEx liver: Spearman rho minus the median rho of nearest protein-coding gene pairs in the same TSS-distance bin (0-5, 5-20, 20-50, 50-100 kb); enhancer-like chromatin at the dominant TSS in HepG2 (highest-TPM isoform for loci, Ensembl canonical for controls; ENCODE fold change over control, mean over TSS +- 1 kb; H3K27ac ENCFF795ONN, H3K4me1 ENCFF576YVM, H3K4me3 ENCFF284FVP): log2((H3K4me1 + 1) / (H3K4me3 + 1)) at active TSSs (H3K27ac >= 2-fold), inactive TSSs ranked lowest; low polysome association (low isoform MRL percentile) |
| M4 Trans-acting RNA, nuclear | Abundance (log2 TPM HepG2, `m4_abundance`: a trans-acting RNA needs copies); HepG2 eCLIP binding by nuclear RBPs (see above) |
| M5 Trans-acting RNA, cytoplasmic | Cytoplasmic enrichment (the negative nuclear log2 ratio); HepG2 eCLIP binding by cytoplasmic RBPs. Loci not expressed in HepG2 and without fractionation data have no M5 data (coverage 0). |
| M6 Small-RNA host | snoRNA, scaRNA or miRNA genes within the locus span on the same strand (GENCODE plus UCSC wgRna, as GENCODE lacks e.g. the GAS5 snoRNAs); snoRNA-host-gene features |

## Step 5. Positive controls and tests

Controls are processed by the same code as the candidates (added as extra loci/ORFs).

| Score | Controls | Expected (`controls.PANEL`) |
|:---|:---|:---|
| D1 HCC | HULC (up), LINC00261 (down) | log2FC_paired in the stated direction and FDR < 0.05; in the replication cohorts combined r in that direction and FDR < 0.05 |
| D1 chronic liver disease | COL1A1 | positive association with fibrosis stage (Hoang r > 0, p < 0.05) |
| D2 | BIRC5, CYTOR | HR of overall survival > 1 |
| D3 | PVT1 (CN-driven, 8q24 with MYC) | rho_CN > 0, FDR < 0.05, and D3 > 0 |
| D4 | PNPLA3, TM6SF2 | genome-wide significant liver-trait GWAS association (D4 >= 7.3, p < 5e-8) |
| D5 | HULC (positive), LINC01090 (negative, symbol collision) | HULC D5 percentile >= 0.9; LINC01090 <= 1 liver paper (no LINC01089 papers credited) |
| D6 | HNF1A-AS1 | liver-specific (D6 >= 0.75) |
| Disease score, fitted disease score, HCC score | HULC, H19 | > 90th percentile of the background |
| Chronic liver disease score | COL1A1, PNPLA3, TM6SF2 | > 90th percentile of the background |
| Actionability | GPC3 144-152 (FVGEFFTDV, HLA-A*02:01 vaccine peptide); GPC3 | presentation percentile <= 2; tumour specificity > 0 |
| M1 | MTLN (mitoregulin, former LINC00116), SPAAR (LINC00961) | top class |
| M2 | GAS5, ZFAS1 | M2 score >= 0.6 (M6 also high) |
| M3 | HNF1A-AS1 (HNF1A), LINC00261 (FOXA2) | top class |
| M4 | MALAT1, NEAT1 | top class |
| M5 | NORAD | top class |
| M6 | SNHG1, MIR22HG | top class |
| QC paralog | an expressed functional paralog (to be fixed from the read data); a known mis-mapping case | `expressed_paralog`; `likely_artefact` |

`controls.PANEL` is the one definition of these checks: `report` evaluates its conditions and writes the
expectation text from them. A check with a missing value (control not quantifiable, not tested, or not a locus of
the run) is "not evaluated", never "fail".

Disease and actionability controls are scored as GENCODE genes (`ctrl:<gene>`, union of GENCODE exons) against
the lincRNA background; mechanism controls are scored through their `lincorfz` locus (best-ranked ORF), so
controls that are not lincRNA loci of the run (MTLN, SPAAR) cannot be evaluated for M1; "top class" is the class
with the highest standardized score z. `controls.tsv` lists every expectation with its value and pass / fail / not
evaluated. The mechanism controls are feature checks: they show that the features of a class are defined as
intended (several pass on annotation features alone, e.g. divergent position or small-RNA host status, and M4 on
the very high abundance of MALAT1 and NEAT1), not that the mechanism itself is detected.

Tests: pytest unit tests on synthetic data for every statistic and for the scoring arithmetic; the control panel
is checked on real data by the `report` step. Controls that fail on real data are reported, not fixed by tuning
weights.

## Weight fitting

Reported alternative to the default weights (`weights` step, part of every run): non-negative weights for D1-D4 and D6 (D5 excluded, as the
labels come from the literature). Labels: established liver-disease lncRNAs = background loci with >= 5
non-retracted liver papers (D5 counts); negatives = background loci with an HGNC-style symbol and no liver paper
(loci named only by an Ensembl, StringTie or NONCODE id cannot have papers and are not labelled); loci in between are not
labelled. Model: logistic regression on the component percentiles (missing = 0.5) with non-negative coefficients and
an L2 penalty chosen by 5-fold stratified cross-validated AUC among 0.01, 0.1, 1 and 10 (a coarse log grid); the refitted coefficients scaled to sum 1 are the
weights. `weights.tsv` reports them with the cross-validated AUC and the AUC of the default weights on the same loci;
`disease_score_fitted` (`disease_score_fitted.tsv`, joined into the ORF tables) is reported next to the default
score, never in its place. Single-component AUCs are reported.

Result on run v2: 164 positives, 998 negatives. Fitted weights D3 0.46, D1 0.36, D4 0.18, D2 0, D6 0;
cross-validated AUC 0.73, against 0.64 for the default weights on the same loci (similar for every penalty).
Single-component AUCs: D1 0.64, D3 0.63, D4 0.55, D2 0.48, D6 0.45. The labels mark lncRNAs that were studied in
liver, which is not the same as being disease-relevant, so the fit is reported as a sensitivity analysis: the
default weights stay, and `disease_score_fitted` is reported, not used.

## Outputs

Per ORF (`lincdisease.orfs.tsv`; candidates only in `lincdisease.candidates.tsv`, ordered by `disease_score_noD5`
and then the lincorfz `priority_score`; ORFs with a CRISPR screen
phenotype (screen value >= 1, see M1), candidates or not, in `lincdisease.screen_supported.tsv` with the screen
columns, sorted by screen value, confidence tier and disease score): QC flags and tier,
copy-attribution confidence, disease components, disease score (+ direction, coverage, no-D5 version),
HCC and chronic liver disease statistics per cohort,
actionability score, mechanism class scores and standardized scores, dominant class. Per locus (`disease.loci.tsv`):
all locus-level statistics, percentiles and components, for loci and controls. Step tables: `qc.orfs.tsv`,
`peptides.orfs.tsv`, `support.orfs.tsv`, `mechanism.orfs.tsv`, `coexpression.loci.tsv`, `encode.loci.tsv`,
`nmd.loci.tsv` (NMD sensitivity, M2), `D3.somatic.tsv`,
`D1.chronic.tsv`, `D1.hcc_replication.tsv`, `D4.gwas.tsv`, `D4.eqtl.tsv`, `D5.literature.tsv`, `homology.orfs.tsv`, `peptide_features.orfs.tsv`, `hla.orfs.tsv`, `hla.control.tsv`, `eclip.loci.tsv`,
`chromatin.loci.tsv`, `weights.tsv`, `disease_score_fitted.tsv`, `quant/` (regions, features, samples, sums, size
factors), `quant_sra/`, `quant_hcc/`, `quant_nmd/` (samples, sums, size factors); control panel: `controls.tsv`.

Identifiers: `feature_id` is a lincorfz locus (gene id) or a control `ctrl:<gene>`; sample ids are
`<source>:<id>` (`tcga:`, `gtex:`, `sra:`, `govaere:<run>`); `orf_id` and `locus_id` come from lincorfz.
Column prefixes: `stat_<X>` the statistic of component X, `pct_<X>` its percentile among the background,
`comp_<X>` the component entering a score (= pct, or a sub-composite for D1 and D6), `raw_<feature>` a mechanism
feature before ranking, `M<k>_score` / `M<k>_coverage` / `M<k>_z` the class score, its coverage and its
standardized score (`length_bin`, `expression_bin` the bins they were ranked in), `*_direction` the sign summary, `*_cell` statistics adjusted for cell-type markers, `<cohort>_r`
/ `_p` / `_n` per-cohort statistics.
Browsable table: `lincdisease.table.html` (all ORFs, coloured and filterable by category: candidate, of interest (config `html.of_interest`), screen, other; sortable; summary scores
first, then D1-D6, M1-M6 and the detail; with an index of the scores). Methods with the derivations of every
metric: `docs/methods.js` (builds `lincdisease_Methods.docx`). Location: `out_dir` of `config/config.yaml`,
relative to the root directory given at run time.

## Running

The steps run as a Snakemake workflow (`workflow/`, started with `run.sh`; usage in `workflow/README.md`) that
starts from the finished lincorfz run and the candidate table. Public data are downloaded once into the data
directory by the `fetch` step on the login node (compute nodes have no network) and recorded in its
`downloads.tsv`; files without a stable URL are listed with their source in `resources.MANUAL` and checked by
`fetch`. Unit tests: `cd lincdisease && <env>/bin/python -m pytest -q`.
