#!/bin/bash
#SBATCH --job-name=gov_locus_cov
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
# usage: govaere_locus_cov.sh GOVAERE_DIR REGIONS_TSV GENOME_FA MODULES
# Coverage of the quantified loci (quant/regions.tsv: feature_id, chrom, start0, end1) in every complete Govaere
# CRAM (cram/RUN.ok) that has no GOVAERE_DIR/locus_cov/RUN.tsv yet: summed per-base depth of primary alignments
# (samtools bedcov, unstranded like the recount3 bigWigs of the other cohorts) per feature. Re-run it as more runs
# finish. The reads of the regions are first extracted in one pass (multi-region iterator): region-by-region
# queries of the archive CRAMs (large slices) take ~0.25 s each. As a job array (sbatch --array=0-N), task i
# takes every (N+1)-th run starting at the i-th.
set -euo pipefail
gd=$1; regions=$2; fa=$3
set +u   # Lmod scripts use unset variables
type module > /dev/null 2>&1 || source "$MODULESHOME/init/bash"
module load $4
set -u
mkdir -p "$gd/locus_cov"
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
bed=$W/regions.bed
awk -F'\t' -v OFS='\t' 'NR>1{print $2, $3, $4, $1}' "$regions" | sort -k1,1 -k2,2n > "$bed"
task=${SLURM_ARRAY_TASK_ID:-0}; n_tasks=${SLURM_ARRAY_TASK_COUNT:-1}; i=-1
for ok in "$gd"/cram/*.ok; do
  i=$((i + 1)); [ $((i % n_tasks)) -eq "$task" ] || continue
  run=$(basename "$ok" .ok)
  out=$gd/locus_cov/$run.tsv
  [ -s "$out" ] && continue
  samtools view -b -M -L "$bed" -T "$fa" -o "$W/sub.bam" "$gd/cram/$run.cram"
  samtools index "$W/sub.bam"
  samtools bedcov -Q 0 -j "$bed" "$W/sub.bam" |
    awk -F'\t' -v OFS='\t' 'BEGIN{print "feature_id", "ncov"} {s[$4] += $5} END{for (f in s) print f, s[f]}' \
    > "$out.tmp.$$"
  mv "$out.tmp.$$" "$out"
  echo "$run done"
done
