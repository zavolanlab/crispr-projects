#!/bin/bash
# usage: govaere_align.sh RUN FQ1 FQ2 STAR_INDEX GENOME_FA OUT_DIR THREADS MODULES
# STAR -> coordinate-sorted CRAM 3.1 archive (all alignments incl. multimappers, unmapped reads kept with their
# mates, soft clips, original qualities) + STAR gene counts. Writes OUT_DIR/cram/RUN.ok only if the CRAM holds
# every input read (primary records = 2 x input pairs); otherwise exits with an error.
set -euo pipefail
run=$1; fq1=$2; fq2=$3; index=$4; fa=$5; out=$6; threads=$7
set +u   # Lmod scripts use unset variables
type module > /dev/null 2>&1 || source "$MODULESHOME/init/bash"
module load $8
set -u
mkdir -p "$out/cram" "$out/counts" "$out/star_logs"
W=${TMPDIR:-/tmp}/govaere_$run; rm -rf "$W"; mkdir -p "$W"
trap 'rm -rf "$W"' EXIT
STAR --runThreadN "$threads" --genomeDir "$index" --readFilesCommand zcat --readFilesIn "$fq1" "$fq2" \
     --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs \
     --outSAMattributes NH HI AS nM NM MD --outSAMattrRGline ID:$run SM:$run \
     --quantMode GeneCounts --limitBAMsortRAM 20000000000 --outFileNamePrefix "$W/"
samtools view -@ "$threads" -C -T "$fa" --output-fmt-option embed_ref=0 --output-fmt-option version=3.1 \
     --output-fmt-option archive -o "$out/cram/$run.cram" "$W/Aligned.sortedByCoord.out.bam"
samtools index "$out/cram/$run.cram"
cp "$W/ReadsPerGene.out.tab" "$out/counts/$run.ReadsPerGene.tsv"
cp "$W/Log.final.out" "$out/star_logs/$run.Log.final.out"
cp "$W/SJ.out.tab" "$out/star_logs/$run.SJ.out.tab"
n_in=$(awk -F'|' '/Number of input reads/{gsub(/[ \t]/,"",$2); print $2}' "$W/Log.final.out")
n_cram=$(samtools view -c -F 0x900 -T "$fa" "$out/cram/$run.cram")
echo "$run input_pairs $n_in primary_records $n_cram"
[ "$n_cram" -eq $((2 * n_in)) ] || { echo "$run INCOMPLETE CRAM"; exit 1; }
echo "$run complete" > "$out/cram/$run.ok"
