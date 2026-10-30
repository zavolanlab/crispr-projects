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
* peptide homology to annotated proteins (DIAMOND against UniProt and GENCODE translations): parent protein,
  identity, coverage;
* flags already computed by `lincorfz`: antisense overlap with coding exons, small-RNA overlap,
  `resolved_multimapper_dominated`, possible 5'/3' extension of a coding gene;
* `tcga_quantifiable`: TCGA and GTEx RNA-seq is unstranded, so for loci antisense to expressed genes only exonic
  positions not overlapping the sense gene are quantified; loci with too few such positions are flagged.

**Copy-attribution confidence** (0-1), for ORFs in duplicated regions or with homology to annotated proteins:

1. copy-specific expression and translation: fraction of P-sites from uniquely mapped reads, in-frame P-sites on
   positions that differ from the paralogs, uniquely mapped RNA-seq coverage and copy-specific splice junctions;
2. copy-specific peptides: identified peptides (proteomics, immunopeptidomics) covering residues that differ from
   the paralogs (1 if present; missing otherwise, as absence of peptides is not evidence against the copy);
3. selection on the copy: 1 - min(Ka/Ks, 1) between the ORF and its parent/paralog ORFs, ORF integrity since
   duplication;
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
* recount3 SRA projects of liver disease cohorts (MASLD/MASH with histological fibrosis and activity scores; HBV/HCV
  liver tissue; additional HCC tumour/normal cohorts for replication), as available.

Values (`ncov`): summed coverage / number of quantifiable positions / the sample's total base coverage
(recount3 `bc_auc.all_reads_all_bases`) x 1e9. Statistics use log2(ncov + pseudocount).

## Step 3. Disease components (per locus)

| | Component | Measure |
|:---|:---|:---|
| D1 | Disease-associated expression | HCC: rank-biserial correlation r of TCGA-LIHC primary tumours vs adjacent normals of the same patients (matched pairs, Wilcoxon) and of tumours vs GTEx liver (Mann-Whitney; AUC = (r+1)/2); statistic = mean of the two r (log2FC and BH-FDR reported); replication in independent HCC cohorts. Chronic liver disease: association with fibrosis stage and MASH activity, and disease vs control liver, in MASLD/MASH and viral-hepatitis cohorts. HCC carries 0.7 and chronic liver disease 0.3 of D1. |
| D2 | Prognosis | Cox models of overall survival and progression-free interval in TCGA-LIHC (TCGA-CDR endpoints; one primary tumour per patient) on standardized expression, adjusted for age and stage (I-II vs III-IV); loci expressed in < 20% of tumours are not tested; statistic = mean Wald z of OS and PFI (HR per SD and BH-FDR reported) |
| D3 | Somatic genome | Copy-number-driven expression in cis (CN-expression correlation, gain/loss frequency); promoter-methylation silencing (450K methylation vs expression). In-ORF mutations are reported but not scored (exome capture barely covers lincRNAs). |
| D4 | Germline genetics | GWAS Catalog associations with whitelisted liver traits (EFO labels of HCC, MASLD/MASH, liver fibrosis, cirrhosis, alcoholic liver disease, viral hepatitis, ...; ALT, AST, GGT) within the locus span +- 50 kb (LD-window proxy); statistic = largest -log10 p (0 without association), reported per trait class (disease / enzyme); GTEx liver eQTLs where the gene is annotated |
| D5 | Prior knowledge | Europe PMC, title or abstract: exact symbol or an HGNC previous symbol / alias (up to 5; only aliases of >= 5 characters that contain a digit and are not another gene's approved symbol, as words like ALIEN or short tokens like LEU2 collect unrelated papers), AND an RNA term (guards against symbols that are words), AND a liver term; counts of all, liver, review and retracted papers; statistic = log2(1 + liver papers - retracted). Exact-phrase matching keeps papers on LINC01089 from being credited to LINC01090. Loci named only by an Ensembl or StringTie id count 0. |
| D6 | Liver context | Equal-weight composite of the percentiles of: median expression in GTEx liver; liver specificity (liver median minus the largest median of the other GTEx tissues); max(log2 TPM HepG2, Huh7) from `lincorfz` |

**Disease score** = weighted mean of component percentiles, shrunk toward 0.5 for missing components.
Default weights: D1 0.25, D2 0.25, D3 0.20, D4 0.10, D5 0.10, D6 0.10. Reported with:

* `disease_direction`: from the signs of the informative signed components (D1 HCC, D1 chronic, D2, D3 with
  percentile >= 0.75; D3 signed so that > 0 means expression raised by copy gain or loss of methylation):
  oncogene-like if all positive, suppressor-like if all negative, inconsistent if both, none if no component
  is informative;
* `disease_coverage`, the confidence tier, and `disease_score_noD5` (without prior knowledge).

**Actionability score** (per ORF, separate from the disease score): equal-weight composite of the percentiles
(among all ORFs) of tumour-specific expression of the locus (median in TCGA-LIHC primary tumours minus the
largest median over GTEx tissues except testis), predicted HLA class I binding (MHCflurry; HepG2/Huh7 HLA types
and common alleles), and the number of lincORF-specific HLA-I immunopeptides (`lincorfz` notebook section 7.3).

## Step 4. Mechanism classes (per ORF, with locus-level features)

Features are ranked against all ORFs of the run; locus-level features are copied to every ORF of the locus.
Each class has an absolute evidence score (mean of its feature percentiles, missing features counting 0.5) and a
relative share (its score divided by the sum over classes). Classes are not mutually exclusive; a class is called
dominant only if its share exceeds the next one by >= 0.05.

| Class | Features |
|:---|:---|
| M1 Encoded peptide | Translation strength (pooled log10 Bayes factor) and fraction of the transcript's P-sites in the ORF (`lincorfz`); identified peptides (proteomics peptides absent from canonical proteins, plus lincORF-specific immunopeptides); independent translation evidence (GENCODE Ribo-seq ORFs, OpenProt); PhyloCSF weighted by copy-attribution confidence (0.5 + (percentile - 0.5) x confidence); ORF vs non-ORF exon conservation; transmembrane, signal-peptide, Pfam features; low predicted disorder; length; hits in published sORF frameshift CRISPR screens |
| M2 Translation-coupled RNA regulation | Stop codon > 50 nt upstream of the last exon junction (NMD substrate); upregulation upon UPF1 knock-down (ENCODE HepG2); number of other translated ORFs (moderate/strong evidence) on the locus |
| M3 Cis-regulation of a neighbour (RNA or DNA element) | Divergent (head-to-head, TSSs within 5 kb) or antisense (exon overlap) position next to a protein-coding gene; co-expression with the neighbour across TCGA and GTEx (vs gene pairs at matched distances); nuclear enrichment (ENCODE HepG2 fractionation); enhancer-like chromatin at the TSS in HepG2; low polysome association (low isoform MRL percentile) |
| M4 Trans-acting RNA, nuclear | Nuclear enrichment; high copy number (HepG2 TPM); HepG2 eCLIP binding by chromatin/transcription-related RBPs; no cis signature |
| M5 Trans-acting RNA, cytoplasmic | Cytoplasmic enrichment; HepG2 eCLIP binding by cytoplasmic RBPs; ribosome association without one dominant ORF (MRL percentile x (1 - largest ORF share of P-sites on the locus)); miRNA-site density (low weight) |
| M6 Small-RNA host | snoRNA, scaRNA or miRNA genes within the locus span on the same strand (GENCODE plus UCSC wgRna, as GENCODE lacks e.g. the GAS5 snoRNAs); snoRNA-host-gene features |

## Step 5. Positive controls and tests

Controls are processed by the same code as the candidates (added as extra loci/ORFs).

| Score | Controls | Expected |
|:---|:---|:---|
| D1 HCC | HULC (up), LINC00261 (down) | correct direction, FDR < 0.05 |
| D1 chronic liver disease | to be fixed after cohort selection (a fibrosis-associated gene, e.g. COL1A1) | positive association with fibrosis stage |
| D2 | BIRC5, CYTOR | HR > 1 |
| D3 | PVT1 (CN-driven, 8q24 with MYC), MEG3 (methylation-silenced) | CN-expression correlation; methylation-expression anticorrelation |
| D4 | PNPLA3, TM6SF2 | liver-trait GWAS association |
| D5 | HULC (positive), LINC01090 (negative, symbol collision) | HULC high; no LINC01089 papers credited to LINC01090 |
| D6 | HNF1A-AS1 | liver-specific |
| Disease score | HULC, H19 | > 90th percentile of the background |
| Actionability | GPC3 144-152 (HLA-A*02:01 vaccine peptide) | predicted binder, tumour-specific |
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

Tests: pytest unit tests on synthetic data for every statistic and for the scoring arithmetic; the control panel
is checked on real data by the `report` step. Controls that fail on real data are reported, not fixed by tuning
weights.

## Weight fitting

Optional, after the default run: non-negative weights for D1-D4 and D6 (D5 excluded to avoid circularity) fitted by
penalized logistic regression of an independent set of established liver-disease lncRNAs against the background,
with cross-validation. Reported as an alternative score with its cross-validated AUC next to the default score.

## Outputs

Per ORF (`lincdisease.orfs.tsv`; candidates only in `lincdisease.candidates.tsv`): QC flags and tier,
copy-attribution confidence, disease components, disease score (+ direction, coverage, no-D5 version),
actionability score, mechanism class scores and relative shares, dominant class. Per locus (`disease.loci.tsv`):
all locus-level statistics, percentiles and components, for loci and controls. Step tables: `qc.orfs.tsv`,
`mechanism.orfs.tsv`, `D4.gwas.tsv`, `D5.literature.tsv`, `quant/`; control panel: `controls.tsv`. Location: `${WF_RUNS_DIR}Liver_lincRNA_ORFs_Claude/disease_mechanism/<run>/`.
