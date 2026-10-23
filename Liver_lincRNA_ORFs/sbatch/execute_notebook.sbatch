#!/bin/bash
#SBATCH --job-name=lincorfz_notebook
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
# --output/--error are given on the sbatch command line (the directory must exist beforehand):
#   sbatch --output=<dir>/notebook_%j.out --error=<dir>/notebook_%j.err \
#          execute_notebook.sbatch <notebook.ipynb> <conda env prefix with nbclient + ipykernel>
#
# Executes a Jupyter notebook in place on a compute node (outputs and figures are stored in it).

set -euo pipefail

NB=$(realpath "$1")
ENV_PREFIX=$(realpath "$2")
cd "$(dirname "${NB}")"

"${ENV_PREFIX}/bin/python" - "${NB}" <<'EOF'
import sys
import nbformat
from nbclient import NotebookClient

path = sys.argv[1]
nb = nbformat.read(path, as_version=4)
client = NotebookClient(nb, timeout=1500, kernel_name="python3",
                        resources={"metadata": {"path": "."}})
try:
    client.execute()
finally:
    nbformat.write(nb, path)
print("executed", path)
EOF
