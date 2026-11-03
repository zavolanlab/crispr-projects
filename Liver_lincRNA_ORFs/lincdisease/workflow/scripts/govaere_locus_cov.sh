#!/bin/bash
#SBATCH --job-name=gov_locus_cov
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
# usage: govaere_locus_cov.sh GOVAERE_DIR REGIONS_TSV GENOME_FA MODULES [RUN]
# Coverage of the quantified loci (quant/regions.tsv: feature_id, chrom, start0, end1) in Govaere CRAMs: summed
# per-base depth of primary alignments (samtools bedcov, unstranded like the recount3 bigWigs of the other cohorts)
# per feature, one row per feature, in GOVAERE_DIR/locus_cov/RUN.tsv. With RUN (workflow rule govaere_locus_cov):
# that run only. Without RUN (by hand, e.g. sbatch --array=0-N): every complete CRAM (cram/RUN.ok) without a
# complete locus_cov file; task i of an array takes every (N+1)-th run starting at the i-th.
# GOVAERE_DIR/locus_cov/regions.tsv keeps the regions the files were computed on; the script stops if REGIONS_TSV
# differs from it (delete locus_cov/ to recompute everything on the new regions).
# The reads of the regions are first extracted in one pass (multi-region iterator): region-by-region queries of
# the archive CRAMs (large slices) take ~0.25 s each.
set -euo pipefail
shopt -s nullglob
gd=$1; regions=$2; fa=$3; modules=$4; only=${5:-}
set +u   # Lmod scripts use unset variables
type module > /dev/null 2>&1 || source "$MODULESHOME/init/bash"
module load $modules
set -u
mkdir -p "$gd/locus_cov"
saved=$gd/locus_cov/regions.tsv
if [ -e "$saved" ]; then
  # features in both tables whose regions differ (as expression.changed_features)
  changed=$(awk -F'\t' 'FNR == 1 {next} NR == FNR {a[$1] = a[$1] " " $2 ":" $3 "-" $4; next}
                        {b[$1] = b[$1] " " $2 ":" $3 "-" $4}
                        END {n = 0; for (f in b) if ((f in a) && a[f] != b[f]) n++; print n}' "$saved" "$regions")
  [ "$changed" -eq 0 ] || { echo "regions of $changed features differ from $saved: delete $gd/locus_cov"; exit 1; }
else
  cut -f1-4 "$regions" > "$saved.tmp.$$" && mv "$saved.tmp.$$" "$saved"
fi
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
bed=$W/regions.bed
awk -F'\t' -v OFS='\t' 'NR>1{print $2, $3, $4, $1}' "$regions" | sort -k1,1 -k2,2n > "$bed"
n_feat=$(cut -f4 "$bed" | sort -u | wc -l)

one_run() {
  local run=$1 out=$gd/locus_cov/$1.tsv
  samtools view -b -M -L "$bed" -T "$fa" -o "$W/sub.bam" "$gd/cram/$run.cram"
  samtools index "$W/sub.bam"
  samtools bedcov -Q 0 -j "$bed" "$W/sub.bam" |
    awk -F'\t' -v OFS='\t' 'BEGIN{print "feature_id", "depth_sum"} {s[$4] += $5} END{for (f in s) print f, s[f]}' \
    > "$out.tmp.$$"
  mv "$out.tmp.$$" "$out"
  echo "$run done"
}

if [ -n "$only" ]; then
  one_run "$only"
  exit 0
fi
task=${SLURM_ARRAY_TASK_ID:-0}; n_tasks=${SLURM_ARRAY_TASK_COUNT:-1}; i=-1
for ok in "$gd"/cram/*.ok; do
  i=$((i + 1)); [ $((i % n_tasks)) -eq "$task" ] || continue
  run=$(basename "$ok" .ok)
  out=$gd/locus_cov/$run.tsv
  [ -s "$out" ] && [ "$(($(wc -l < "$out") - 1))" -eq "$n_feat" ] && continue
  one_run "$run"
done
