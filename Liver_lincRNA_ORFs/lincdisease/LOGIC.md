# lincdisease: disease relevance and functional mechanisms of translated lincRNA ORFs

This file describes the current logic of the `lincdisease` module. It is kept as one clean, up-to-date
specification: when the logic changes, this file is edited in place (no change log; history is in git).

## Scope and inputs

* **Input**: a table of ORFs produced by `lincorfz` (same columns as `results/final/lincORFs.ranked.tsv`), taken as
  is. The current input is the manually curated `lincORFs.subselection.ranked.tsv.csv` of run `v2` (46 ORFs, 33 loci).
* **Background**: all expressed lincRNA loci and their ORFs from the same `lincorfz` run. Every measure is also
  computed for the background, and scores are expressed as percentiles against it (0.5 = a typical expressed lincRNA).
* **Disease focus**: liver disease, with hepatocellular carcinoma (HCC) as the main but not the only focus
  (MASLD/MASH and fibrosis, viral hepatitis, liver-function traits).
* **Two levels**: disease evidence attaches mostly to the lincRNA **locus**; mechanism evidence separates what the
  **RNA** does from what the **peptide** does. Locus and ORF scores are reported side by side and never multiplied.

## General rules

1. Each component is a statistic converted to an empirical percentile among the background (signed statistics
   are ranked on magnitude; the sign is kept separately as a direction).
2. Missing data never counts as zero: a score is shrunk toward 0.5 in proportion to the weight of its missing
   components, and is reported with its `coverage` (fraction of component weight with data).
3. Quality checks produce flags and a confidence tier; they do not remove ORFs and do not change the scores.
4. Default weights are fixed in advance. Fitted weights (see "Weight fitting") are reported as an alternative,
   never silently substituted.

## Step 1. Quality-check layer (per ORF)

Flags:

* overlap with segmental duplications (UCSC `genomicSuperDups`, with identity `fracMatch`) and with GRCh38 regions
  that are falsely duplicated relative to T2T-CHM13;
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

1. copy-specific expression and translation: fraction of P-sites from uniquely mapped reads, in-frame P-sites on
   positions that differ from the paralogs, uniquely mapped RNA-seq coverage and copy-specific splice junctions;
2. copy-specific peptides: identified peptides (proteomics, immunopeptidomics) covering residues that differ from
   the paralogs (1 if present; missing otherwise, as absence of peptides is not evidence against the copy).
   Implemented as a peptide absent from canonical proteins and contaminants that maps to this lincORF only;
3. selection on the copy: 1 - min(Ka/Ks, 1), Nei-Gojobori Ka/Ks (Jukes-Cantor corrected) between the ORF CDS and
   the CDS of the best GENCODE hit, on the codon alignment implied by the protein alignment; missing with < 20 aligned
   codons, dS = 0 (identical copies) or saturated p-distances. Paralogous lincORFs without an annotated parent are not
   compared;
4. alignment attribution: age of the duplication, (0.98 - identity) / 0.08 clipped to 0-1 (identity = segmental
   duplication `fracMatch`, else protein identity to the homolog). For old duplications the multi-species
   alignment likely contains true orthologs of the copy; for very recent ones (> 98% identity) other species
   probably carry the parent.

The confidence is the mean of the available parts (1 for ORFs outside duplications and without homologs).

Conservation is **not** discarded in duplicated regions: if the copy is expressed and translated, conservation
indicates evolutionary support for its function. PhyloCSF enters the peptide score weighted by the
copy-attribution confidence.

**Confidence tier**: `clean`; `caution` (antisense overlap with coding exons, small-RNA overlap,
multimapper-dominated, possible coding-gene extension, not quantifiable in TCGA/GTEx); `expressed_paralog`
(duplicated or homologous, with copy-specific evidence); `likely_artefact` (duplicated or homologous, and the
translation evidence cannot be attributed to the copy: < 50% unique P-sites or multimapper-dominated, no
copy-specific peptide, and no selection on the copy).

## Step 2. Quantification in patient and tissue data

Many loci are de novo assemblies without a GENCODE gene, so gene-level public tables cannot be used. Coverage is
read from recount3 base-level coverage bigWigs, per locus (disease evidence attaches to the locus), over its
**quantifiable positions**: the union of the exons of its isoforms minus the exons of every GENCODE gene on
either strand that is not the locus itself (public RNA-seq is unstranded). Controls are measured the same way on
the union of their GENCODE exons. Loci with < 100 quantifiable positions are flagged `tcga_quantifiable = False`.
Samples:

* TCGA-LIHC (all tumours and adjacent normals),
* GTEx (all normal liver samples; 10 random samples of every other tissue),
* recount3 SRA liver disease cohorts (`quant_sra/`): Hoang 2019 (GSE130970, 78 NAFLD biopsies with fibrosis stage and
  NAFLD activity score), Suppli 2019 (GSE126848, 57: normal-weight, obese, NAFL, NASH), Boldanova 2017 (GSE84346,
  untreated biopsies, HCV vs non-HCV). Govaere 2020 (GSE135251) is not in recount3. Phenotypes from the recount3
  SRA `sample_attributes`.

Values (`ncov`): summed coverage / number of quantifiable positions / the sample's total base coverage
(recount3 `bc_auc.all_reads_all_bases`) x 1e9. Statistics use log2(ncov + pseudocount) (chronic disease cohorts: rescaled by gene-based size factors, see D1).

## Step 3. Disease components (per locus)

| | Component | Measure |
|:---|:---|:---|
| D1 | Disease-associated expression | HCC: rank-biserial correlation r of TCGA-LIHC primary tumours vs adjacent normals of the same patients (matched pairs, Wilcoxon) and of tumours vs GTEx liver (Mann-Whitney; AUC = (r+1)/2); statistic = mean of the two r (log2FC and BH-FDR reported); replication in independent HCC cohorts. Chronic liver disease: ncov rescaled by median-of-ratios size factors of the recount3 gene sums (expressed protein-coding genes, top 1% by abundance left out), because AUC scaling lets the hepatocyte transcripts that dominate liver RNA set the scale and their falling share with injury lifted 74% of loci with fibrosis; per cohort, partial rank correlation with the phenotype (ranks residualized on log input reads, % uniquely mapped reads and exonic fraction; t test with n - 2 - k df): Hoang vs fibrosis stage (vs NAFLD activity score reported), Suppli NAFL + NASH vs normal-weight + obese, Boldanova HCV vs non-HCV; a feature detected (ncov > 0) in < 20% of a cohort's samples is not tested there; statistic = weighted mean of the scored cohort statistics (fibrosis 2/3, HCV 1/3), combined p by Stouffer (signed z, weights sqrt(n)), BH-FDR. Suppli is reported, not scored: its cases differ technically from its controls (fewer uniquely mapped reads, lower exonic fraction; the 15 samples with the lowest exonic fraction are all cases), and after rescaling and adjustment 40% of loci keep p < 0.01 (median r +0.14; Hoang 9%, 0.00). Parent genes: a locus whose expression follows a parent gene (homology step; median Spearman over cohorts >= 0.5) that is itself associated with chronic disease (FDR < 0.05, same sign) is `parent_driven`, and its cohort statistics are recomputed with the parent's expression as an extra covariate (processed pseudogenes pick up reads of parents such as KRT8, KRT18, YWHAZ that rise in injured liver). Reported, not scored: `D1_chronic_cell`, the same statistic additionally adjusted for liver cell-type marker scores (hepatocyte, stellate, cholangiocyte, macrophage, lymphocyte, endothelial; MacParland et al. 2018 markers, COL1A1, KRT8, KRT18 left out), i.e. the part not explained by a change in cell mix. HCC carries 0.7 and chronic liver disease 0.3 of D1. |
| D2 | Prognosis | Cox models of overall survival and progression-free interval in TCGA-LIHC (TCGA-CDR endpoints; one primary tumour per patient) on standardized expression, adjusted for age and stage (I-II vs III-IV); loci expressed in < 20% of tumours are not tested; statistic = mean Wald z of OS and PFI (HR per SD and BH-FDR reported) |
| D3 | Somatic genome | TCGA-LIHC primary tumours (GDC GRCh38 data from the UCSC Xena GDC hub). Copy number: masked CN segment value (log2 ratio) at the midpoint of the locus span; Spearman rho of CN vs expression; gain / loss frequency (> +0.3 / < -0.3). Methylation: 450K probes within -1500..+500 nt of the TSS of the dominant isoform (loci: highest TPM in HepG2/Huh7; GENCODE controls: Ensembl canonical transcript; windows around every isoform TSS would pull gene-body probes of long multi-isoform genes into the promoter), mean beta per sample; Spearman rho of beta vs expression; delta_beta = median beta in adjacent normals - median in tumours. Signed parts: d3_cn = max(rho_CN, 0) x (gain - loss frequency); d3_meth = max(-rho_meth, 0) x clip(delta_beta / 0.2, -1, 1); statistic = the part with the larger magnitude. Both rho are partial Spearman correlations adjusted for tumour purity (ABSOLUTE, TCGA PanCanAtlas; rank residuals after linear regression on purity rank, tumours with a purity call), as non-tumour cells shift expression, CN amplitude and methylation together; unadjusted rho are reported. In-ORF mutations are reported but not scored (exome capture barely covers lincRNAs). |
| D4 | Germline genetics | GWAS Catalog associations with whitelisted liver traits (EFO labels of HCC, MASLD/MASH, liver fibrosis, cirrhosis, alcoholic liver disease, viral hepatitis, ...; ALT, AST, GGT) within the locus span +- 50 kb (LD-window proxy); statistic = largest -log10 p (0 without association), reported per trait class (disease / enzyme). GTEx v8 liver eQTLs for loci with GENCODE genes (v26 ids without version): eGene q < 0.05, lead variant, and a link to a whitelisted GWAS hit if a significant eQTL variant is the hit or the lead eQTL lies within 10 kb (distance proxy, no LD panel); reported (`D4.eqtl.tsv`, `D4_eqtl_supported`), not scored, as only annotated loci can have GTEx eQTLs |
| D5 | Prior knowledge | Europe PMC, title or abstract: exact symbol or an HGNC previous symbol / alias (up to 5; only aliases of >= 5 characters that contain a digit and are not another gene's approved symbol, as words like ALIEN or short tokens like LEU2 collect unrelated papers), AND an RNA term (guards against symbols that are words), AND a liver term; counts of all, liver, review and retracted papers; statistic = log2(1 + liver papers - retracted). Exact-phrase matching keeps papers on LINC01089 from being credited to LINC01090. Loci named only by an Ensembl or StringTie id count 0. |
| D6 | Liver context | Equal-weight composite of the percentiles of: median expression in GTEx liver; liver specificity (liver median minus the largest median of the other GTEx tissues); max(log2 TPM HepG2, Huh7) from `lincorfz` |

**Disease score** = weighted mean of component percentiles, shrunk toward 0.5 for missing components.
Default weights: D1 0.25, D2 0.25, D3 0.20, D4 0.10, D5 0.10, D6 0.10. Reported with:

* `disease_direction`: from the signs of the informative signed components (D1 HCC, D1 chronic, D2, D3 with
  percentile >= 0.75; D3 signed so that > 0 means expression raised by copy gain or loss of methylation):
  oncogene-like if all positive, suppressor-like if all negative, inconsistent if both, none if no component
  is informative;
* `disease_coverage`, the confidence tier, and `disease_score_noD5` (without prior knowledge).

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
(among all ORFs) of tumour-specific expression of the locus (median in TCGA-LIHC primary tumours minus the
largest median over GTEx tissues except testis), predicted HLA class I presentation (MHCflurry 2 presentation
percentile <= 2, distinct 8-11-mer binders over HepG2 A*02:01, A*24:02, B*35:14, B*51:01, C*04:01, C*16:02, Huh7
A*11:01, B*54:01, C*01:02 (Cellosaurus) and 16 common alleles, each scored on its own), and the number of lincORF-specific HLA-I immunopeptides (`lincorfz` notebook section 7.3).

## Step 4. Mechanism classes (per ORF, with locus-level features)

Features are ranked against all ORFs of the run; locus-level features are copied to every ORF of the locus.
ENCODE HepG2 features use stranded data, so they are measured on the exons of the feature minus the exons of other
GENCODE genes on the same strand only (`quant/regions.tsv` also removes antisense exons, which unstranded RNA-seq
needs but stranded data do not). Fractionation uses stranded signal bigWigs of uniquely mapped reads (two
replicates, averaged), as mean signal per position on the locus strand: nuclear enrichment = log2 (nucleus + 0.01) / (cytosol +
0.01) from CeFra-seq, cytoplasmic enrichment its negative; missing if both conditions have < 0.05 signal per base.
HepG2 eCLIP (M4, M5): reproducible (IDR) peaks of 105 RBPs (ENCODE), log2 fold enrichment >= 3 and -log10 p >= 3,
on the locus strand. RBP groups by HepG2 immunofluorescence (Van Nostrand et al. 2020, Suppl. Data 1): nuclear =
only nuclear compartments (17 RBPs), cytoplasmic = only cytoplasmic compartments (15); RBPs seen in both (61) or not
imaged (12) count only in the total. Statistic per group: log2[((n_group + f) / (n_bound + 1)) / f], f = the group's
share of the 105 RBPs (over-representation among the bound RBPs, independent of expression and length); missing if
no RBP binds.
Each class has an absolute evidence score (mean of its feature percentiles, missing features counting 0.5) and a
relative share (its score divided by the sum over classes). Classes are not mutually exclusive; a class is called
dominant only if its share exceeds the next one by >= 0.05.

| Class | Features |
|:---|:---|
| M1 Encoded peptide | Translation strength (pooled log10 Bayes factor) and fraction of the transcript's P-sites in the ORF (`lincorfz`); identified peptides (proteomics peptides absent from canonical proteins, plus lincORF-specific immunopeptides); independent translation evidence: same stop codon as an ORF of the GENCODE Ribo-seq ORF catalog (phase 2 comprehensive, v45), +1 if it matches a phase 1 ORF with peptide support (OpenProt not used); PhyloCSF weighted by copy-attribution confidence (0.5 + (percentile - 0.5) x confidence); ORF vs non-ORF exon conservation: mean PhyloP100way over the ORF minus the mean over the exonic positions of its transcript outside the ORF, both restricted to the locus' quantifiable positions (Step 2: exons of other GENCODE genes on either strand removed, as PhyloP is strand-blind and antisense coding exons or miRNA hairpins would lend their conservation to the ORF; missing with < 30 ORF or < 30 other positions), weighted by copy-attribution confidence like PhyloCSF; length; sequence features = 0.5 x [Kyte-Doolittle TM segment (19-aa window mean >= 1.6), signal-peptide-like N-terminus (heuristic: hydrophobic 8-mer, mean >= 2.0 without D/E/K/R, starting at residues 2-20, small residues at -3/-1 of a cleavage site <= 40; open predictors only, SignalP/TMHMM are licensed), or Pfam-A domain (hmmsearch, gathering thresholds)] + 0.5 x (1 - fraction of residues disordered, metapredict >= 0.5). Published CRISPR knock-out screens, matched by the same stop codon: (a) Prensner et al. 2021 (553 non-canonical ORFs; primary screen in 8 cell lines including HepG2, secondary tiling screen of the hits): 2 = viability phenotype specific to the ORF (secondary screen: specific to the ORF, or to the ORF and one other transcript region), 1 = phenotype not resolved by the secondary screen, 0 = no phenotype or a phenotype attributed to the transcript, the genomic region or an adjacent coding gene; the HepG2 median sgRNA log2 fold change is reported. (b) Chen et al. 2020 (~2,300 ORFs, K562 and iPSC; coordinates lifted from hg19 with the UCSC chain): 1 = hit in either cell line, 0 = screened without a hit. Their hit calls are not in the supplementary tables, and their stated null (non-targeting sgRNAs) makes nearly every ORF significant, because Cas9 cutting alone costs a phenotype of about -0.015 (their estimate); hits are re-called as BH-FDR < 0.05 in a one-sided Mann-Whitney test of the ORF's sgRNA phenotypes against the non-targeting sgRNAs shifted by -0.015, and a top-3 sgRNA phenotype <= -0.02 (808 hits in either cell line, 146 in both; the paper reports 401 in both, so the re-call is stricter). Feature: Prensner's secondary-screen verdict where it exists (the only test of ORF specificity), otherwise the larger of the two values; missing if neither screen tested the ORF. |
| M2 Translation-coupled RNA regulation | Stop codon > 50 nt upstream of the last exon junction (NMD substrate); upregulation upon NMD inhibition (missing: the ENCODE HepG2 UPF1 shRNA knock-down lowers UPF1 mRNA by only ~30% and does not raise canonical NMD targets such as GAS5, ATF4, GADD45B, TBL2, so it is not used); number of other translated ORFs (moderate/strong evidence) on the locus |
| M3 Cis-regulation of a neighbour (RNA or DNA element) | Divergent (head-to-head, TSSs within 5 kb) or antisense (exon overlap) position next to a protein-coding gene; co-expression with the nearest protein-coding gene (TSS distance <= 100 kb, either strand; recount3 gene sums) across TCGA-LIHC and GTEx liver: Spearman rho minus the median rho of nearest protein-coding gene pairs in the same TSS-distance bin (0-5, 5-20, 20-50, 50-100 kb); nuclear enrichment (ENCODE HepG2 fractionation); enhancer-like chromatin at the dominant TSS in HepG2 (highest-TPM isoform for loci, Ensembl canonical for controls; ENCODE fold change over control, mean over TSS +- 1 kb; H3K27ac ENCFF795ONN, H3K4me1 ENCFF576YVM, H3K4me3 ENCFF284FVP): log2((H3K4me1 + 1) / (H3K4me3 + 1)) at active TSSs (H3K27ac >= 2-fold), inactive TSSs ranked lowest; low polysome association (low isoform MRL percentile) |
| M4 Trans-acting RNA, nuclear | Nuclear enrichment; high copy number (HepG2 TPM); HepG2 eCLIP binding by nuclear RBPs (see above); no cis signature |
| M5 Trans-acting RNA, cytoplasmic | Cytoplasmic enrichment; HepG2 eCLIP binding by cytoplasmic RBPs; ribosome association without one dominant ORF (MRL percentile x (1 - largest ORF share of P-sites on the locus)); miRNA-site density (low weight) |
| M6 Small-RNA host | snoRNA, scaRNA or miRNA genes within the locus span on the same strand (GENCODE plus UCSC wgRna, as GENCODE lacks e.g. the GAS5 snoRNAs); snoRNA-host-gene features |

## Step 5. Positive controls and tests

Controls are processed by the same code as the candidates (added as extra loci/ORFs).

| Score | Controls | Expected |
|:---|:---|:---|
| D1 HCC | HULC (up), LINC00261 (down) | correct direction, FDR < 0.05 |
| D1 chronic liver disease | COL1A1 | positive association with fibrosis stage (Hoang rho > 0, p < 0.05) |
| D2 | BIRC5, CYTOR | HR > 1 |
| D3 | PVT1 (CN-driven, 8q24 with MYC), MEG3 (reported methylation-silenced), GSTP1 (promoter hypermethylation in HCC; added post hoc, see below) | PVT1: rho_CN > 0 (FDR < 0.05) and D3 > 0; MEG3, GSTP1: rho_meth < 0 (FDR < 0.05) |
| D4 | PNPLA3, TM6SF2 | liver-trait GWAS association |
| D5 | HULC (positive), LINC01090 (negative, symbol collision) | HULC high; no LINC01089 papers credited to LINC01090 |
| D6 | HNF1A-AS1 | liver-specific |
| Disease score | HULC, H19 | > 90th percentile of the background |
| HCC score | HULC, H19 | > 90th percentile of the background |
| Chronic liver disease score | COL1A1, PNPLA3, TM6SF2 | > 90th percentile of the background |
| Actionability | GPC3 144-152 (FVGEFFTDV, HLA-A*02:01 vaccine peptide) | presentation percentile <= 2; GPC3 tumour-specific |
| M1 | MTLN (mitoregulin, former LINC00116), SPAAR (LINC00961) | top class |
| M2 | GAS5, ZFAS1 | high (M6 also high) |
| M3 | HNF1A-AS1 (HNF1A), LINC00261 (FOXA2) | top class |
| M4 | MALAT1, NEAT1 | top class |
| M5 | NORAD | top class |
| M6 | SNHG1, MIR22HG | top class |
| QC paralog | an expressed functional paralog (to be fixed from the read data); a known mis-mapping case | `expressed_paralog`; `likely_artefact` |

Disease and actionability controls are scored as GENCODE genes (`ctrl:<gene>`, union of GENCODE exons) against
the lincRNA background; mechanism controls are scored through their `lincorfz` locus (best-ranked ORF), so
controls that are not lincRNA loci of the run (MTLN, SPAAR) cannot be evaluated for M1. `controls.tsv` lists every
expectation with its value and pass / fail / not evaluated.

MEG3 is imprinted and, in TCGA-LIHC tumours, positively correlated with the methylation of its canonical promoter (rho 0.37 after purity adjustment), so its expectation, derived from cell lines and small cohorts, does not hold in this cohort; it stays in the panel as a reported failure. GSTP1 was chosen after this was seen, as the only one of five classic HCC methylation markers (CDKN2A, RASSF1, APC, SFRP1, GSTP1) anticorrelated with expression; it shows that the statistic detects methylation silencing, not a pre-registered expectation.

Tests: pytest unit tests on synthetic data for every statistic and for the scoring arithmetic; the control panel
is checked on real data by the `report` step. Controls that fail on real data are reported, not fixed by tuning
weights.

## Weight fitting

Optional, after the default run (`weights` step): non-negative weights for D1-D4 and D6 (D5 excluded, as the
labels come from the literature). Labels: established liver-disease lncRNAs = background loci with >= 5
non-retracted liver papers (D5 counts); negatives = background loci with an HGNC-style symbol and no liver paper
(loci named only by an Ensembl or StringTie id cannot have papers and are not labelled); loci in between are not
labelled. Model: logistic regression on the component percentiles (missing = 0.5) with non-negative coefficients and
an L2 penalty chosen by 5-fold stratified cross-validated AUC; the refitted coefficients scaled to sum 1 are the
weights. `weights.tsv` reports them with the cross-validated AUC and the AUC of the default weights on the same loci;
`disease_score_fitted` (`disease_score_fitted.tsv`, joined into the ORF tables) is reported next to the default
score, never in its place. Single-component AUCs are reported.

Result on run v2 (184 positives, 1,111 negatives): the fit puts 0.82 on D3 and 0.18 on D4 (cross-validated AUC
0.75; default weights 0.59); without D3 the AUC is 0.53. D3 separates the labels by ascertainment, not disease
biology: established lncRNAs are older annotations covered by 450K probes (78% of positives have promoter probes,
31% of negatives; the probe count alone gives AUC 0.81), and their CN-expression correlation is stronger. Literature
labels therefore cannot fit the weights here; `disease_score_fitted` is reported but not used, and the default
weights stay.

## Outputs

Per ORF (`lincdisease.orfs.tsv`; candidates only in `lincdisease.candidates.tsv`; ORFs with a CRISPR screen
phenotype (screen value >= 1, see M1), candidates or not, in `lincdisease.screen_supported.tsv` with the screen
columns, sorted by screen value, confidence tier and disease score): QC flags and tier,
copy-attribution confidence, disease components, disease score (+ direction, coverage, no-D5 version),
HCC and chronic liver disease statistics per cohort,
actionability score, mechanism class scores and relative shares, dominant class. Per locus (`disease.loci.tsv`):
all locus-level statistics, percentiles and components, for loci and controls. Step tables: `qc.orfs.tsv`,
`peptides.orfs.tsv`, `support.orfs.tsv`, `mechanism.orfs.tsv`, `coexpression.loci.tsv`, `encode.loci.tsv`, `D3.somatic.tsv`,
`D1.chronic.tsv`, `D4.gwas.tsv`, `D4.eqtl.tsv`, `D5.literature.tsv`, `homology.orfs.tsv`, `peptide_features.orfs.tsv`, `hla.orfs.tsv`, `hla.control.tsv`, `eclip.loci.tsv`,
`chromatin.loci.tsv`, `weights.tsv`, `disease_score_fitted.tsv`, `quant/`, `quant_sra/`; control panel: `controls.tsv`.
Browsable table: `lincdisease.table.html` (candidates and screen-supported ORFs; sortable; summary scores
first, then D1-D6, M1-M6 and the detail; with an index of the scores). Methods with the derivations of every
metric: `docs/methods.js` (builds `lincdisease_Methods.docx`). Location: `out_dir` of `config/config.yaml`,
relative to the root directory given at run time.

## Running

The steps run as a Snakemake workflow (`workflow/`, started with `run.sh`; usage in `workflow/README.md`) that
starts from the finished lincorfz run and the candidate table. Public data are downloaded once into the data
directory by the `fetch` step on the login node (compute nodes have no network); files without a stable URL are
listed with their source in `resources.MANUAL` and checked by `fetch`.
