#!/bin/bash
#SBATCH --job-name=lincorfz_driver
#SBATCH --qos=1day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# --output/--error are given on the sbatch command line (the directory must exist beforehand,
# otherwise SLURM fails silently):
#   sbatch --output=<rundir>/slurm_logs/driver_%j.out --error=<rundir>/slurm_logs/driver_%j.err \
#          run_lincorfz.sbatch <params.json> <rundir> <path/to/lincorfz>
#
# Runs the Nextflow driver of lincorfz inside a SLURM job; every pipeline task is in turn
# submitted by Nextflow as its own SLURM job (executor 'slurm'). Compute nodes have no internet
# access, so the conda environments must have been pre-built with lincorfz/setup_conda_envs.sh.

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "usage: sbatch [...] run_lincorfz.sbatch <params.json> <rundir> <lincorfz_dir>" >&2
    exit 1
fi
PARAMS=$(realpath "$1")
RUNDIR=$(realpath -m "$2")
PIPELINE_DIR=$(realpath "$3")

set +u
source ~/.bashrc
conda activate nextflow
set -u

export NXF_OFFLINE=true
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_ANSI_LOG=false

mkdir -p "${RUNDIR}"
cd "${RUNDIR}"   # .nextflow/ cache and .nextflow.log live here, which makes -resume work

nextflow run "${PIPELINE_DIR}/main.nf" \
    -params-file "${PARAMS}" \
    -profile conda \
    -resume
