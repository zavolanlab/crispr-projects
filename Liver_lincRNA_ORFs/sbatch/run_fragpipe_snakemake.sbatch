#!/bin/bash
#SBATCH --job-name=fragpipe_driver
#SBATCH --qos=1day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# --output/--error are given on the sbatch command line (the directory must exist beforehand):
#   sbatch --output=<proteomics_dir>/slurm_logs/driver_%j.out --error=<proteomics_dir>/slurm_logs/driver_%j.err \
#          run_fragpipe_snakemake.sbatch <fragpipe-snakemake dir> <target file> [snakemake conda env]
#
# Runs the Snakemake driver of zavolanlab/fragpipe-snakemake (configured via its config.yaml)
# inside a SLURM job; every rule is submitted as its own SLURM job (profiles/slurm). The rule
# conda environments must have been created beforehand on a node with internet access:
#   snakemake --use-conda --conda-create-envs-only --cores 1 <target file>
set -euo pipefail
if [ $# -lt 2 ]; then
    echo "usage: sbatch [...] run_fragpipe_snakemake.sbatch <fragpipe-snakemake dir> <target file> [conda env]" >&2
    exit 1
fi
WORKFLOW_DIR=$(realpath "$1")
TARGET=$(realpath -m "$2")
ENV_NAME=${3:-snakemake}
set +u
source ~/.bashrc
conda activate "${ENV_NAME}"
set -u
export PYTHONNOUSERSITE=1
cd "${WORKFLOW_DIR}"
snakemake --profile profiles/slurm "${TARGET}"
