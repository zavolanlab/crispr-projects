# lincdisease workflow

Annotates the lincorfz ORFs with disease relevance and mechanism classes (see `../LOGIC.md`). It starts from a
finished lincorfz run and does not rerun any part of it. Everything is driven by `../run.sh`, so you do not need
to know Snakemake to use it.

## Root directory and inputs
All paths in `../config/config.yaml` are relative to a **root directory**, given at run time with
`run.sh --root DIR` or taken from the current directory. No paths are stored in git. Expected layout:

```
<root>/
  inputs/lincorfz_results     -> lincorfz results/ (ranked ORFs, lincRNA transcripts, final/ candidate table)
  inputs/lincorfz_peptides    -> directory with proteomics/ and peptidomics/ from the lincorfz notebook
  inputs/PhyloP100way.bw      -> PhyloP100way bigWig
  conda_envs/                 lincdisease, lincdisease_peptide (from ../envs/), snakemake
  lincdisease_data/           download cache (about 40 GB; filled by the fetch step)
  disease_mechanism/v2/       outputs; the workflow runs here
  external_data/govaere2020/  Govaere FASTQ (temporary), CRAMs, gene counts
```

The three `inputs/` entries are symbolic links to wherever the lincorfz run lives, for example:

```
mkdir -p inputs
ln -s <lincorfz run>/results inputs/lincorfz_results
ln -s <lincorfz project dir> inputs/lincorfz_peptides
ln -s <PhyloP100way.bw> inputs/PhyloP100way.bw
```

To start from another lincorfz run or candidate table, change the links (or the entries in `config.yaml`) and
use a new `out_dir`.

## Daily use
Run these from the root directory on a login node (`L` = the `lincdisease/` directory of the repository):

```
L/run.sh plan            # what would run; changes nothing. L/run.sh plan -v also gives the reason per step
L/run.sh start           # run in the background; you can log out
L/run.sh status          # running or not, progress lines and errors, Slurm jobs
L/run.sh stop            # stop and cancel the Slurm jobs
L/run.sh redo disease    # rerun a step (rule name) and everything that depends on it
L/run.sh start govaere   # Govaere 2020 FASTQ -> CRAM + gene counts
```

From elsewhere, add `--root DIR` before the command. A step reruns by itself when one of its input tables
changes, or when one of the lincdisease modules it uses is edited. Changes in `__main__.py` and new downloads do
not trigger reruns; use `run.sh redo STEP` for those.

Final outputs in the output directory: `lincdisease.table.html` (candidates and screen-supported ORFs, sortable,
with an index of the scores), `lincdisease.candidates.tsv`, `lincdisease.orfs.tsv`,
`lincdisease.screen_supported.tsv`, `controls.tsv`.

## Where to look when something fails
- `run.sh status` shows the error lines. The full log is `<out_dir>/logs/workflow/latest.log`; the failing step's
  own log is printed there too (`show-failed-logs`).
- The log of each step is `<out_dir>/logs/<step>.log`; the Slurm log is `<out_dir>/logs/slurm/<step>/<jobid>.log`.
- "Directory cannot be locked": an earlier run was killed (e.g. the login session ended uncleanly). Check with
  `run.sh status` that nothing is running, then `run.sh unlock` and start again.
- A Slurm job killed for memory is retried once with twice the memory. If it fails again, raise its `mem(...)`
  value in `Snakefile`.
- A failed step does not stop independent steps (`keep-going`); after fixing the cause, `run.sh start` resumes
  and reruns only what is missing or failed.

## How it runs
- `Snakefile` lists one rule per `python -m lincdisease <step>` with its input and output tables;
  `rules/govaere.smk` holds the Govaere alignment.
- Compute nodes have no network. Steps that download (`fetch`, `manifest`, `chronic_manifest`, `quant`,
  `quant_sra`, `literature`, `govaere_fastq`) run on the login node inside the Snakemake process. All other steps
  are submitted with plain `sbatch` (`profile/config.yaml`), and their state is read with `sacct`
  (`scripts/slurm_status.sh`).
- `fetch` downloads every public resource with a stable URL into `data_dir`. It also checks the files that were
  placed by hand (`resources.MANUAL`: journal supplements, GENCODE Ribo-seq ORFs, Pfam, GTEx eQTLs, MHCflurry
  models, proteomes) and lists each missing file with its source.
- The recount3 streaming (`quant`, `quant_sra`) appends one row per sample. An interrupted run keeps what it
  has, and a rerun fetches only the missing samples. Done markers are in `<out_dir>/flags/`.
- Govaere: the FASTQ files are temporary. Snakemake deletes them once their CRAM passes the completeness check
  (`cram/<run>.ok`). A download waits while `max_fastq_gb` of FASTQ is on disk.

## Software
Snakemake 9.27.1 with the cluster-generic executor 1.0.9, in its own environment
(`<root>/conda_envs/snakemake`; set `SNAKEMAKE_ENV` to use another one).
Rebuild it with:

```
module load Miniforge3
mamba create -p <prefix> -c conda-forge -c bioconda python=3.12 snakemake-minimal=9.27.1 \
      snakemake-executor-plugin-cluster-generic=1.0.9
```
