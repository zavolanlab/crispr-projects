#!/bin/bash
# usage: govaere_fetch.sh RUNS_TSV RUN FASTQ_DIR MAX_GB [MAX_WAIT_H]
# Download the two ENA FASTQ files of RUN (login node) and check their md5. Waits first while FASTQ_DIR holds
# MAX_GB or more, so that downloads do not run far ahead of the alignments; gives up after MAX_WAIT_H hours
# (default 6: the FASTQ of failed alignments is not deleted and would otherwise block all downloads).
set -euo pipefail
runs_tsv=$1; run=$2; dir=$3; max_gb=$4; max_wait_h=${5:-6}
mkdir -p "$dir"
waited=0
while [ "$(du -s --block-size=1G "$dir" | cut -f1)" -ge "$max_gb" ]; do
  if [ "$waited" -ge $((max_wait_h * 3600)) ]; then
    echo "$(date) $dir still holds >= $max_gb GB of FASTQ after ${max_wait_h} h; check failed alignments"; exit 1
  fi
  echo "$(date) $dir holds >= $max_gb GB of FASTQ; waiting"; sleep 120; waited=$((waited + 120))
done
line=$(awk -F'\t' -v r="$run" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["run_accession"]==r{print $h["fastq_ftp"]"\t"$h["fastq_md5"]}' "$runs_tsv")
[ -n "$line" ] || { echo "$run not in $runs_tsv"; exit 1; }
IFS=$'\t' read -r urls md5s <<< "$line"
IFS=';' read -ra U <<< "$urls"; IFS=';' read -ra M <<< "$md5s"
for i in "${!U[@]}"; do
  f=$dir/$(basename "${U[$i]}")
  ok=0
  for attempt in 1 2 3; do
    wget -q -c -O "$f" "https://${U[$i]}" || true
    if echo "${M[$i]}  $f" | md5sum -c --quiet; then ok=1; break; fi
    echo "$(date) $f: md5 mismatch (attempt $attempt), downloading again"; rm -f "$f"; sleep 60
  done
  [ $ok = 1 ] || { echo "$f: download failed"; exit 1; }
  echo "$(date) ok $f"
done
