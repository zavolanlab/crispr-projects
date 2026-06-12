# Ribosome profiling reveals extensive translational alterations in steatohepatitis-associated hepatocellular carcinoma

This repository contains the jupyter notebook and auxiliary scripts related to the analysis of evolutionary conservation and mutations associated with novel ORFs (nORFs) identified in Ribo-seq data from tumor and non-tumor hepatocellular carcinoma (HCC) samples.

On sciCORE HPC, running jupyter notebook on a computational node is nicely enabled by [OnDemand service](https://docs.scicore.unibas.ch/HPC%20Cluster/interactivecomputing/#open-ondemand-ood).

**Jupyter Notebook** is used for interactive analysis and visualization, and for automatic generation of SLURM scripts that should be run on the HPC cluster.

## Repository Structure

```text
.
├── CRISPR_Projects-current.ipynb       # Main master Jupyter notebook for downstream analysis
├── CRISPR_projects.template.env        # Template for required environment variables/paths
├── install/
│   └── environment.yaml                # Conda environment specification for the Jupyter notebook
└── WF/                                 # Snakemake Workflow Engine
    ├── Snakefile-prepare-faster        # Pipeline Step 1: RNA-seq data prep (alignment, etc.)
    ├── Snakefile-quantification-faster # Pipeline Step 2: Quantification (gene counts, PAQR quantification of APA usage)
    ├── config.template.yaml            # Template configuration for Snakemake parameters
    ├── envs/                           # Conda environments isolated for specific Snakemake rules
    ├── profile/                        # SLURM execution profile for the HPC
    └── scripts/                        # Python and R scripts utilized by both Snakemake and Jupyter
```

## Quick Start & Setup

To ensure strict reproducibility and security, this project uses `.env` files to manage all absolute paths (data directories, genome annotations, etc.). **Do not hardcode paths into the Python or Snakemake files.**

### 1. Clone the Repository
Clone this repository into your local user space (`$HOME`):
```bash
git clone https://github.com/zavolanlab/crispr-projects.git
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

### 3. Install the conda environment with zavolab_pyutils
Analysis in the notebook is largely based on the functions from [zavolab_pyutils](https://github.com/zavolanlab/zavolab_pyutils/tree/dev) repository.
Follow the instruction from that repo "Developer Setup from source, with conda environment".
Use the created conda environment "zavolab_pyutils" to execute the Jupyter Notebook.

## Executing the Workflows

The heavy lifting is divided into two separate Snakemake workflows located in the `WF/` directory.

Configuration of the workflows (i.e. creation of input .tsv with sample specification and .yaml config is done **inside** the jupyter notebook)

`Bash` commands are also prepared inside the jupyter notebook. They should be further copied into command line and executed.

**On an HPC cluster like sciCORE**, workflows should be executed on a **login** node. Snakemake further automatically submits jobs to computational nodes.

## Downstream Analysis
Once the Snakemake workflows are complete, all results are routed to the shared group directories defined in your `.env` file. 

Use respective sections of `CRISPR_Projects-current.ipynb` to analyze the outputs. 

The notebook automatically loads your `.env` paths using `python-dotenv`, allowing it to dynamically locate all workflow results, figures, and metadata regardless of where you cloned this repository.