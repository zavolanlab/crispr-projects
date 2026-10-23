#!/bin/bash
# Pre-build the conda environments used by lincorfz.
#
# sciCORE compute nodes have no internet access, so environments cannot be created
# by Nextflow inside SLURM jobs. Run this script ONCE on a login node (it only downloads
# and installs packages); afterwards the pipeline activates the pre-built environments
# from <conda_envs_dir>, which must match the `conda_envs_dir` parameter of the run.
#
# The Python environment additionally gets Sanity_utils (private repository nimwegenLab/Sanity_utils;
# SSH access to GitHub is needed, or set SANITY_UTILS_SRC to a local clone) and the C++ Sanity
# binary, built with a portable architecture flag so that it runs on all cluster nodes.
#
# usage: bash setup_conda_envs.sh <conda_envs_dir>
#        SANITY_UTILS_SRC=~/libs/Sanity_utils bash setup_conda_envs.sh <conda_envs_dir>
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: bash $0 <conda_envs_dir>" >&2
    exit 1
fi
ENVS_DIR=$(realpath -m "$1")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$ENVS_DIR"

# mamba is much faster at solving; fall back to conda
SOLVER=$(command -v mamba || command -v conda)

for yml in "$SCRIPT_DIR"/envs/*.yml; do
    name=$(basename "$yml" .yml)
    prefix="$ENVS_DIR/$name"
    if [ -d "$prefix/conda-meta" ]; then
        echo "[skip] $name already exists at $prefix"
        continue
    fi
    echo "[create] $name -> $prefix"
    "$SOLVER" env create -y -p "$prefix" -f "$yml"
done
PY="$ENVS_DIR/lincorfz_py"
if [ ! -x "$PY/bin/Sanity" ]; then
    SRC=${SANITY_UTILS_SRC:-"git+ssh://git@github.com/nimwegenLab/Sanity_utils.git"}
    echo "[install] Sanity_utils from $SRC and the Sanity binary -> $PY"
    PYTHONNOUSERSITE=1 "$PY/bin/pip" install --no-deps "$SRC"
    CONDA_PREFIX="$PY" PATH="$PY/bin:$PATH" "$PY/bin/sanity-utils-install-sanity" \
        --prefix "$PY" --arch-flag=-march=x86-64-v2
fi
echo "done: environments in $ENVS_DIR"
