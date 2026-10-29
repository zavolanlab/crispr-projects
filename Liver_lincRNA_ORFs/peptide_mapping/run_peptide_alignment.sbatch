#!/bin/bash
#SBATCH --job-name=peptide_align_driver
#SBATCH --qos=1day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# --output/--error are given on the sbatch command line (the directory must exist beforehand):
#   sbatch --output=<rundir>/slurm_logs/driver_%j.out --error=<rundir>/slurm_logs/driver_%j.err \
#          run_peptide_alignment.sbatch <NF_pre-rRNA dir> <config> <rundir> [nextflow conda env]
#
# Runs the Nextflow driver of the RNA_oxidation alignment workflow (NF_pre-rRNA,
# splicing_trimmed_alignment_all_samples.nf: adapter trimming, STAR, minimum-NM filtering and
# multimapper resolution) on the peptide-encoding RNA reads; every process is its own SLURM job.
# The config is used exclusively (nextflow -C), as in the original launcher of that workflow.
set -euo pipefail
if [ $# -lt 3 ]; then
    echo "usage: sbatch [...] run_peptide_alignment.sbatch <NF_pre-rRNA dir> <config> <rundir> [conda env]" >&2
    exit 1
fi
NF_DIR=$(realpath "$1")
CONFIG=$(realpath "$2")
RUNDIR=$(realpath -m "$3")
ENV_NAME=${4:-nextflow}
set +u
source ~/.bashrc
conda activate "${ENV_NAME}"
set -u
export NXF_OFFLINE=true
export NXF_OPTS="-Xms1g -Xmx3g"
export NXF_ANSI_LOG=false
mkdir -p "${RUNDIR}"
cd "${RUNDIR}"   # .nextflow/ cache and .nextflow.log live here, which makes -resume work
nextflow -C "${CONFIG}" run "${NF_DIR}/splicing_trimmed_alignment_all_samples.nf" -resume
