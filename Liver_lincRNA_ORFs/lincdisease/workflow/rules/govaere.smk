"""
Govaere et al. 2020 NAFLD cohort (GSE135251, PRJNA558102): ENA FASTQ -> STAR -> CRAM 3.1 archive + STAR gene
counts -> coverage of the quantified loci (locus_cov/<run>.tsv), which the `chronic` rule reads for the runs in
GOVAERE_USED (all runs of the run table except config govaere.exclude). The FASTQ is a temporary file: Snakemake
deletes it once its CRAM is complete, and a download waits while max_fastq_gb of FASTQ is on disk. Target
`govaere`: all CRAMs, without the downstream steps.
"""
import pandas as pd

GOV = config["govaere"]
GOV_DIR, GOV_REF = GOV["dir"], GOV["reference"]
GOV_RUNS_TSV = str(PKG / GOV["runs"])
GOV_RUNS = pd.read_csv(GOV_RUNS_TSV, sep="\t")["run_accession"].tolist()
GOVAERE_EXCLUDE = GOV.get("exclude") or []        # a list: written space-separated into the shell command
GOVAERE_USED = [r for r in GOV_RUNS if r not in GOVAERE_EXCLUDE]
SCRIPTS = str(Path(workflow.basedir) / "scripts")
GENOME = f"{GOV_REF}/GRCh38.primary_assembly.genome.fa"

wildcard_constraints:
    run=r"[SED]RR\d+",


rule govaere:
    """All CRAMs (target `govaere`), without the downstream steps."""
    input: expand(f"{GOV_DIR}/cram/{{run}}.ok", run=GOV_RUNS)


rule govaere_fastq:
    output: temp(f"{GOV_DIR}/fastq/{{run}}_1.fastq.gz"), temp(f"{GOV_DIR}/fastq/{{run}}_2.fastq.gz")
    resources: govaere_download=1
    log: f"{GOV_DIR}/logs/fetch/{{run}}.log"
    shell: "bash {SCRIPTS}/govaere_fetch.sh {GOV_RUNS_TSV} {wildcards.run} {GOV_DIR}/fastq {GOV[max_fastq_gb]} > {log} 2>&1"


rule star_index:
    """STAR index of the GENCODE v42 primary assembly (2x76 nt reads: overhang 75); the FASTA is the CRAM reference."""
    input: gz=ancient(f"{DATA}/homology/GRCh38.primary_assembly.genome.fa.gz"),
           gtf=ancient(f"{DATA}/gencode.v42.annotation.gtf.gz")
    output: fa=GENOME, fai=GENOME + ".fai", sa=f"{GOV_REF}/star_v42_oh75/SA"
    threads: 8
    resources: mem_mb=64000, runtime=240
    log: f"{GOV_DIR}/logs/star_index.log"
    shell:
        """
        (set +u
         type module > /dev/null 2>&1 || source "${{MODULESHOME}}/init/bash"
         module load {GOV[modules]}
         set -u
         zcat {input.gz} > {output.fa}
         samtools faidx {output.fa}
         zcat {input.gtf} > $TMPDIR/annotation.gtf
         STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {GOV_REF}/star_v42_oh75 \
              --genomeFastaFiles {output.fa} --sjdbGTFfile $TMPDIR/annotation.gtf --sjdbOverhang 75 \
              --outTmpDir $TMPDIR/star_idx_tmp) > {log} 2>&1
        """


rule govaere_align:
    input: fq1=f"{GOV_DIR}/fastq/{{run}}_1.fastq.gz", fq2=f"{GOV_DIR}/fastq/{{run}}_2.fastq.gz",
           index=ancient(f"{GOV_REF}/star_v42_oh75/SA"), fa=ancient(GENOME), fai=ancient(GENOME + ".fai")
    output: cram=f"{GOV_DIR}/cram/{{run}}.cram", crai=f"{GOV_DIR}/cram/{{run}}.cram.crai",
            counts=f"{GOV_DIR}/counts/{{run}}.ReadsPerGene.tsv", star_log=f"{GOV_DIR}/star_logs/{{run}}.Log.final.out",
            sj=f"{GOV_DIR}/star_logs/{{run}}.SJ.out.tab", ok=f"{GOV_DIR}/cram/{{run}}.ok"
    threads: 8
    resources: mem_mb=64000, runtime=240     # 48 GB ran out of memory on some runs
    priority: 10          # align what is on disk before downloading more
    log: f"{GOV_DIR}/logs/align/{{run}}.log"
    shell:
        "bash {SCRIPTS}/govaere_align.sh {wildcards.run} {input.fq1} {input.fq2} {GOV_REF}/star_v42_oh75 {input.fa} "
        "{GOV_DIR} {threads} '{GOV[modules]}' > {log} 2>&1"


rule govaere_locus_cov:
    """Coverage of the quantified loci (quant/regions.tsv) in one CRAM. The regions are an `ancient` input, as for
    the recount3 sums: locus_cov/regions.tsv keeps the regions the files were computed on, and the script and the
    chronic step stop if they have changed (delete locus_cov/ to recompute)."""
    input: ok=f"{GOV_DIR}/cram/{{run}}.ok", regions=ancient("quant/regions.tsv"), fa=ancient(GENOME)
    output: f"{GOV_DIR}/locus_cov/{{run}}.tsv"
    threads: 2
    resources: mem_mb=mem(4), runtime=30
    log: f"{GOV_DIR}/logs/locus_cov/{{run}}.log"
    shell: "bash {SCRIPTS}/govaere_locus_cov.sh {GOV_DIR} {input.regions} {input.fa} '{GOV[modules]}' {wildcards.run} "
           "> {log} 2>&1"
