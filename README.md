# CRISPR Projects Analysis and Pipelines

This repository contains the computational workflows and downstream analysis notebooks related to CRISPR projects in Zavolan lab. 
It utilizes a hybrid approach: **Snakemake** for robust, scalable data processing on HPC clusters (sciCORE), and **Jupyter Notebooks** for interactive downstream analysis and visualization.

# Current state
Currently, we use RNA-seq and Ribo-seq data to identify targets for CRISPR perturbations
We make UCSC trackhub to interactively visualize these data using [UCSC Genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A6534517-6538371)

## Repository Structure

```text
.
├── CRISPR_Projects-current.ipynb       # Main master Jupyter notebook for downstream analysis
├── CRISPR_projects.template.env        # Template for required environment variables/paths
├── install/
│   └── environment.yaml                # Conda environment specification for the Jupyter notebook
└── WF/                                 # Snakemake Workflow Engine
    ├── Snakefile-prepare-faster        # Pipeline Step 1: RNA-seq data prep (alignment, etc.)
    ├── Snakefile-quantification-faster # Pipeline Step 2: Quantification (gene counts, PAQR APA)
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
git clone [https://github.com/mirono0000/crispr-projects.git](https://github.com/mirono0000/crispr-projects.git)
cd crispr-projects
```

### 2. Configure Environment Paths
You must map the project to your local HPC paths. 
* **If you are a group member on sciCORE:** Simply symlink the master shared configuration into your local repository directory:
  ```bash
  ln -s /scicore/home/zavolan/GROUP/CRISPR_projects/configs/CRISPR_projects.scicore.env .env
  ```
* **If you are an external collaborator:** Copy the template, rename it, and fill in your absolute paths:
  ```bash
  cp CRISPR_projects.template.env .env
  # Open .env and edit the "Base Directories" section to match your system
  ```
*(Note: `*.env` files are ignored by git to protect private cluster paths).*

### 3. Install the Conda Environment
Create and activate the master environment required to run the Jupyter notebook and standard data science libraries (Pandas, UMAP, SciPy, BioPython, etc.):
```bash
conda env create -f install/environment.yaml
conda activate crispr_projects

# Link the environment to your Jupyter kernels
python -m ipykernel install --user --name crispr_projects --display-name "Python (crispr_projects)"
```

---

## Executing the Workflows

The heavy lifting is divided into two separate Snakemake workflows located in the `WF/` directory.

### Step 1: Data Preparation
Processes raw sequencing data, handles stranded alignments, and generates genome-browser-ready tracks (e.g., BigWigs).
```bash
cd WF
# Copy and edit the template config for your specific run
cp config.template.yaml config.yaml 

# Execute via SLURM profile
snakemake -s Snakefile-prepare-faster --profile profile/
```

### Step 2: Quantification
Takes the prepared alignments and quantifies gene expression, transcript usage, and PAQR APA (Alternative Polyadenylation) counts.
```bash
snakemake -s Snakefile-quantification-faster --profile profile/
```

---

## Downstream Analysis
Once the Snakemake workflows are complete, all results are routed to the shared group directories defined in your `.env` file. 

Open `CRISPR_Projects-current.ipynb` to begin analyzing the outputs. The notebook automatically loads your `.env` paths using `python-dotenv`, allowing it to dynamically locate all workflow results, figures, and metadata regardless of where you cloned this repository.