# Liver lincRNA ORFs: automated discovery and prioritization of translated lincRNA ORFs

This sub-project automates the search for translated ORFs in long intergenic non-coding RNAs (lincRNAs)
of hepatic cell lines. It replaces the manual, IGV-based pre-selection that led to the ZFAS1 and FRG1BP
targets of the HCC Ribo-seq paper (siRNA knock-down of the lincRNA, CRISPR pgRNA excision of the ORF).

The analysis consists of

* **`lincorfz/`**: a Nextflow (DSL2) pipeline that does all heavy computations as SLURM jobs on sciCORE,
* **`Liver_lincRNA_ORFs.ipynb`**: the Jupyter notebook that prepares the inputs, writes the parameter
  files, prints the `sbatch` commands and inspects/visualizes the results,
* **`sbatch/run_lincorfz.sbatch`**: the SLURM job that runs the Nextflow driver.

All results are written to `/scicore/home/zavolan/GROUP/Liver_projects/wf_runs/Liver_lincRNA_ORFs_Claude/`.

## What the pipeline produces

The main deliverable is a GTF in hg38 coordinates for IGV, `results/final/lincORFs.igv.gtf`, containing

* every expressed, de novo assembled lincRNA isoform (`transcript_type "lincRNA_isoform"`), and
* every ORF with Ribo-seq support (`transcript_type "lincRNA_ORF"`), drawn on its isoform: exons, `CDS`,
  `start_codon` and `stop_codon` rows. The ORF `transcript_id` is `ORF<rank>_<locus name>_<length>aa`, so
  the rank is visible in IGV; all evidence (P-sites, periodicity p-value, MRL, expression, PhyloCSF,
  amino-acid sequence) is stored as GTF attributes and shown when clicking on the feature.

It comes with an IGV session (`results/final/lincorfz.igv_session.xml`) that loads the GTF next to the
supporting tracks, a ranked table of all ORFs and a protein FASTA.

| file (under `results/`) | content |
|:---|:---|
| `final/lincORFs.igv.gtf` | lincRNA isoforms and ORFs with Ribo-seq evidence (IGV-ready, sorted) |
| `final/lincORFs.ranked.tsv` | all ORFs with >= 1 P-site, all metrics, sorted by evidence class and priority score |
| `final/lincORFs.loci.ranked.tsv` | best ORF per lincRNA locus |
| `final/lincORFs.orfs.bed` | BED12 of the ORFs (colour = Ribo-seq evidence class, score = priority) |
| `final/lincORFs.faa` | ORF protein sequences (e.g. for peptidomics searches) |
| `final/lincorfz.igv_session.xml` | IGV desktop session with all tracks below |
| `polysome_mrl/polysome.mrl.{plus,minus}.bw` | per-nucleotide mean ribosome load, 80S..decasome |
| `polysome_mrl/polysome.total.{plus,minus}.bw` | polysomal coverage (sum of fraction means) |
| `polysome_mrl/polysome.frac{k}.{plus,minus}.bw` | replicate-mean coverage of the fraction with k ribosomes |
| `polysome_mrl/polysome.mrl_rep{r}.{plus,minus}.bw` | per-replicate MRL (replicates with all fractions) |
| `riboseq/psite_bigwigs/pooled_cpm.{plus,minus}.bw` | HepG2 Ribo-seq P-sites, 3 replicates pooled, CPM |
| `riboseq/psite_bigwigs/<sample>.{plus,minus}.bw` | P-site counts per replicate |
| `riboseq/calibration/` | P-site offsets per read class and start/stop metagene plots |
| `lincRNAs/lincRNA_transcripts.{gtf,tsv}` | selected expressed lincRNA isoforms and their annotation |
| `lincRNAs/lincRNA_candidates.all.tsv` | all non-coding assembled transcripts with the reason for (non-)selection |
| `expression/lincRNA_loci.expression.tsv` | pySanity log2 CPM (per sample, per cell line, errors) |
| `polysome_scores/`, `conservation/`, `orfs/` | per-module tables |
| `pipeline_info/` | Nextflow report, timeline, trace, DAG |

## Data

| layer | data (sciCORE) |
|:---|:---|
| bulk RNA-seq | HepG2 and Huh7, 3 replicates each, TruSeq stranded PE (F2R1), PMID 40885748; STAR genome alignments and `SJ.out.tab` in `wf_runs/Asiers_HepG2_Huh7_RNAseq_v1/` |
| Ribo-seq | HepG2 CHX-only, 3 replicates, genomic STAR alignments with probabilistically resolved multimappers (`RNA_oxidation/wf_runs/Juan_pre_rRNA/.../resolution/*.resolved.bam`) |
| polysome profiling | HepG2 fractions 80S, disome ... decasome, spike-in normalized stranded coverage (`wf_runs/polysome_seq_hepg2_inayat_v1/polysome_normalised_biwig_data/`) |
| annotation | GENCODE v42 (comprehensive), NONCODE v6 lncRNAs, HCC-paper lincRNA ORFs (`tables/HCC_ORFs/annotations/nORFs_CDSs.gtf`) |
| conservation | PhyloCSF frame tracks and power, PhyloP 100-/30-way (`external_data/conservation_and_mutations/`) |

## Method

**1. Mean ribosome load (MRL).** For every position *x* and strand, *c<sub>k</sub>(x)* is the mean over replicates
of the spike-in normalized coverage of the fraction carrying *k* ribosomes (80S: *k* = 1 ... decasome: *k* = 10;
the free, 40S, 60S and "heavysome" fractions are not used). Replicate means use the replicates available for
each fraction (the decasome has two). The MRL is *MRL(x) = Σ k·c<sub>k</sub>(x) / Σ c<sub>k</sub>(x)* and is written only
where the summed coverage is at least `mrl_min_total_coverage` (default 10, i.e. roughly 10-20 reads;
one read adds 0.4-1.3 normalized units). The input bigwigs use Ensembl chromosome names; outputs use `chr` names.
For isoforms and ORFs, the MRL is computed from the coverage summed over all exonic positions
(*Σ k·C<sub>k</sub> / Σ C<sub>k</sub>*), and reported as a percentile among Ensembl-canonical protein-coding mRNAs.
Because polysome profiling sediments whole RNA molecules, the isoform-level MRL is the primary score.

**2. Ribo-seq P-sites from genomic alignments.** Reads carry ~3 untemplated nucleotides at the 5' end
(5' soft clip of 3 nt in ~55% of reads) and 3' poly-A remnants. P-site offsets are therefore calibrated
per read class (aligned length x 5' soft-clip length) and counted from the aligned 5' end: for reads around
annotated start codons of canonical CDSs, the offset is the distance with the strongest start-codon peak
that also places the P-sites in the dominant CDS frame; classes with < 45% in-frame P-sites are discarded.
Spliced reads are kept only if all their introns are annotated or supported by >= 2 unique RNA-seq junction
reads (this removes spurious long-gap alignments of short reads).

**3. De novo lincRNA annotation.** StringTie (`--rf`, guided by GENCODE) assembles each RNA-seq sample;
assemblies are merged and re-quantified per sample. An assembled isoform is kept as a lincRNA if it
(same strand) does not overlap exons of CDS-containing transcripts, exons of any transcript of a
protein-coding gene, or introns of protein-coding genes (likely pre-mRNA); has >= 200 nt; has a mean
TPM >= 1 in HepG2 or Huh7 (>= 5 for mono-exonic isoforms); every intron is supported by >= 3 unique
junction reads (STAR `SJ.out.tab`, all samples); and makes up >= 10% of its locus' expression.
Overlaps with GENCODE lncRNAs (e.g. ZFAS1), pseudogenes (e.g. FRG1BP), small RNAs, NONCODE genes,
antisense coding genes and possible UTR extensions of neighbouring coding genes are reported, not filtered.

**4. Expression.** featureCounts (reverse-stranded pairs) on GENCODE genes plus the de novo lincRNA loci,
normalized jointly with pySanity (`apply_sanity_normalization_full_bayesian`, zavolab_pyutils); log2 CPM
per sample and per cell line with errors.

**5. ORFs.** All ATG-initiated ORFs (>= 8 codons) of the selected isoforms are enumerated; for ORFs sharing a
stop codon the most upstream start with in-frame P-site support is chosen. Periodicity is tested as in the
HCC Ribo-seq paper (one-sided Wilcoxon rank-sum of per-codon in-frame vs. each off-frame P-site counts,
combined with Fisher's method), per replicate and pooled. The *island score* is the log2 ratio of the P-site
density in the ORF to that in the rest of the isoform. ORFs found in several isoforms are collapsed.

**6. Conservation.** PhyloCSF is quantified exactly as in `HCC_RiboSeq_paper.ipynb`: each ORF segment is
scored on the PhyloCSF frame track matching its codon phase (cumulative ORF length upstream of the segment;
minus-strand frames counted from the chromosome end), segments overlapping annotated CDSs are removed
(`bedtools subtract`), and the median over all positions is reported (likewise PhyloCSF power, PhyloP).

**7. Ranking.** Ribo-seq evidence classes (pooled data):

| class | periodicity p | in-frame P-sites | in-frame fraction | replicates with p < 0.05 |
|:---|:---|:---|:---|:---|
| strong | < 1e-3 | >= 10 | >= 0.5 | >= 2 |
| moderate | < 1e-2 | >= 5 | >= 0.4 | - |
| weak | < 0.05 | >= 3 | - | - |
| multimapping | any of the above, but < 10% of the ORF's P-sites come from uniquely mapped reads | | | |

The *multimapping* class captures ORFs whose Ribo-seq signal consists of multimapping reads that were assigned
probabilistically; this is typical for processed pseudogenes (e.g. ribosomal protein pseudogenes), whose
reads stem from the parent gene, and also applies to FRG1BP (paralogous FRG1 family). The class based on all
reads and the periodicity of uniquely mapped reads alone are reported in separate columns
(`ribo_evidence_all_reads`, `ribo_evidence_unique_reads`).

Within the classes, ORFs are ordered by a priority score: the weighted mean of percentile ranks of isoform
MRL (0.25), periodicity -log10 p (0.20), island score (0.10), HepG2 expression (0.20) and PhyloCSF (0.10),
the fraction of P-sites from uniquely mapped reads (0.15), plus peptide evidence (0.05) if peptide tables are
given. All thresholds and weights are parameters.

**8. Optional peptide evidence.** With `peptide_tables` (TSV/CSV with a peptide column, e.g. FragPipe
`peptide.tsv` from non-tryptic small-peptidome searches), identified peptides absent from the GENCODE
proteome are matched to the ORF sequences (I/L-equivalent). `lincORFs.faa` can be used as a search database.

## How to run

1. **Conda environments** (once, on a login node, because compute nodes have no internet access):
   ```bash
   bash lincorfz/setup_conda_envs.sh /scicore/home/zavolan/GROUP/Liver_projects/wf_runs/Liver_lincRNA_ORFs_Claude/conda_envs
   ```
   A `nextflow` conda environment with Nextflow >= 24 is also required (`conda create -n nextflow -c bioconda nextflow`).
2. **Inputs and parameters**: run the section "Prepare inputs" of `Liver_lincRNA_ORFs.ipynb`; it writes the sample
   sheets (`inputs/`) and the parameter files (`params/*.json`).
3. **Run** (the notebook prints this command; execute it on a login node):
   ```bash
   RUN=/scicore/home/zavolan/GROUP/Liver_projects/wf_runs/Liver_lincRNA_ORFs_Claude/runs/v1
   mkdir -p $RUN/slurm_logs
   sbatch --output=$RUN/slurm_logs/driver_%j.out --error=$RUN/slurm_logs/driver_%j.err \
       sbatch/run_lincorfz.sbatch <params/v1.json> $RUN $(pwd)/lincorfz
   ```
   Re-submitting the same command resumes the run (`-resume`). Note that Nextflow does not notice edits of
   the scripts in `lincorfz/bin/`: after changing a script, delete the run's `work/` directory (or the work
   directories of the affected tasks) to recompute. A test on one chromosome is configured by
   `"chromosomes": "chr20"` (ZFAS1 and FRG1BP are both on chr20; it takes ~15 min).
4. **Inspect** the results with the "Results" sections of the notebook and in IGV with
   `results/final/lincorfz.igv_session.xml`.

The pipeline can also be run directly (`nextflow run lincorfz/main.nf -params-file params.json -profile conda`);
see `lincorfz/nextflow.config` for all parameters and their defaults.

## Limitations

* Ribo-seq and polysome profiling are available for HepG2 only; Huh7 contributes expression only.
* Multimapping Ribo-seq reads were assigned probabilistically upstream of this pipeline. ORFs whose signal is
  dominated by such reads are put in the `multimapping` class; they are listed in the tables but, by default,
  not in the GTF (set `gtf_min_evidence` to `multimapping` to include them).
* ORFs overlapping snoRNAs or other small RNAs are flagged (`orf_overlaps_small_RNA`): small-RNA fragments
  can contaminate Ribo-seq libraries.
* Mono-exonic assembled transcripts are kept only when highly expressed, but can still be pre-mRNA or
  enhancer RNA fragments; check `n_exons` and `possible_*_extension_of`.
