# Ribosome profiling reveals extensive translational alterations in steatohepatitis-associated hepatocellular carcinoma

This repository contains the jupyter notebook and auxiliary scripts related to the analysis of evolutionary conservation and mutations associated with novel ORFs (nORFs) identified in Ribo-seq data from tumor and non-tumor hepatocellular carcinoma (HCC) samples.

Historically, this analysis is related to CRISPR projects of the Zavolan Lab, hence it's stored in respective repository.

> **Additional project: Liver lincRNA ORFs.** Automated de novo annotation of expressed lincRNAs in hepatic cell lines and prioritization of their translated ORFs (Ribo-seq, polysome profiling, RNA-seq, PhyloCSF), with a Nextflow pipeline and a dedicated notebook. See [Liver_lincRNA_ORFs/README.md](Liver_lincRNA_ORFs/README.md).

On sciCORE HPC, running jupyter notebook on a computational node is nicely enabled by [OnDemand service](https://docs.scicore.unibas.ch/HPC%20Cluster/interactivecomputing/#open-ondemand-ood).

**Jupyter Notebook** is used for interactive analysis and visualization, and for automatic generation of SLURM scripts that should be run on the HPC cluster.

## Repository Structure

```text
.
├── HCC_RiboSeq_paper.ipynb             # Main master Jupyter notebook with analysis
├── CRISPR_projects.template.env        # Template for required environment variables/paths
└── analysis_scripts/                    
    ├── compute_track_metrics.py        # python script for computing various metrics stored in external .bb or .bigwig files
```

## Quick Start & Setup

To ensure strict reproducibility **and security**, this project uses `.env` files to manage all absolute paths (data directories, genome annotations, etc.). **Do not hardcode absolute paths into the code.**

### 1. Clone the Repository
Clone **HCC_Riboseq_paper branch**:
```bash
git clone --branch HCC_Riboseq_paper https://github.com/zavolanlab/crispr-projects.git
cd crispr-projects
```

### 2. Configure Environment Paths
You must map the project to your local HPC paths. 
**First**, copy the template, rename it, and fill in your absolute paths, for example like that:
  ```bash
  cp CRISPR_projects.template.env CRISPR_projects.scicore.env
  # Open .env and edit the "Base Directories" section to match your system
  ```
* **Recommended if you are a group member on sciCORE:** move the `CRISPR_projects.scicore.env` to Project GROUP folder and symlink into your local repository directory:
  ```bash
  ln -s <a file with specified sciCORE paths> CRISPR_projects.scicore.env
  ```
This way `CRISPR_projects.scicore.env` will be automatically accessible by group members but will not be tracked by git.
*(Note: `*.env` files are ignored by git to protect private cluster paths, except the `CRISPR_projects.template.env` file).*

The notebook automatically loads your `.env` paths using `python-dotenv`, allowing it to dynamically locate the files in the system regardless of where you cloned this repository.

### 3. Install the conda environment with zavolab_pyutils
Analysis in the notebook is largely based on the functions from [zavolab_pyutils](https://github.com/zavolanlab/zavolab_pyutils/tree/dev) repository.
Follow the instruction from that repo "Developer Setup from source, with conda environment".
Use the created conda environment "zavolab_pyutils" to execute the Jupyter Notebook.

## Submitting the SLURM jobs and computation-heavy bash commands

`Bash` commands are also prepared inside the jupyter notebook. They should be further copied into command line and executed (either at the login node or at the computational node, see comments in the notebook).

**On an HPC cluster like sciCORE**, commands for downloading public data from web should be executed on a **login** node. Computation-heavy jobs should be executed strictly from computational nodes.