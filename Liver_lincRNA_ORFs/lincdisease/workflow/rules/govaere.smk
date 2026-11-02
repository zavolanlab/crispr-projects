"""
Govaere et al. 2020 NAFLD cohort (GSE135251, PRJNA558102): ENA FASTQ -> STAR -> CRAM 3.1 archive + STAR gene
counts. The FASTQ is a temporary file: Snakemake deletes it once its CRAM is complete, and a download waits while
max_fastq_gb of FASTQ is on disk. Target: `govaere` (not part of `all` until the cohort enters D1 chronic).
"""
import pandas as pd

GOV = config["govaere"]
GD, REF = GOV["dir"], GOV["reference"]
GOV_RUNS_TSV = str(PKG / GOV["runs"])
GOV_RUNS = pd.read_csv(GOV_RUNS_TSV, sep="\t")["run_accession"].tolist()
SCRIPTS = str(Path(workflow.basedir) / "scripts")
GENOME = f"{REF}/GRCh38.primary_assembly.genome.fa"

wildcard_constraints:
    run=r"[SED]RR\d+",
    w=r"\d+",



rule govaere:
    input: expand(f"{GD}/cram/{{run}}.ok", run=GOV_RUNS)


rule govaere_fastq:
    output: temp(f"{GD}/fastq/{{run}}_1.fastq.gz"), temp(f"{GD}/fastq/{{run}}_2.fastq.gz")
    resources: govaere_download=1
    log: f"{GD}/logs/fetch/{{run}}.log"
    shell: "bash {SCRIPTS}/govaere_fetch.sh {GOV_RUNS_TSV} {wildcards.run} {GD}/fastq {GOV[max_fastq_gb]} > {log} 2>&1"


rule star_index:
    """STAR index of the GENCODE v42 primary assembly (2x76 nt reads: overhang 75); the FASTA is the CRAM reference."""
    input: gz=ancient(f"{D}/homology/GRCh38.primary_assembly.genome.fa.gz"),
           gtf=ancient(f"{D}/gencode.v42.annotation.gtf.gz")
    output: fa=GENOME, fai=GENOME + ".fai", sa=f"{REF}/star_v42_oh75/SA"
    threads: 8
    resources: mem_mb=64000, runtime=240
    log: f"{GD}/logs/star_index.log"
    shell:
        """
        (set +u
         type module > /dev/null 2>&1 || source "${{MODULESHOME}}/init/bash"
         module load {GOV[modules]}
         set -u
         zcat {input.gz} > {output.fa}
         samtools faidx {output.fa}
         zcat {input.gtf} > $TMPDIR/annotation.gtf
         STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {REF}/star_v42_oh75 \
              --genomeFastaFiles {output.fa} --sjdbGTFfile $TMPDIR/annotation.gtf --sjdbOverhang 75 \
              --outTmpDir $TMPDIR/star_idx_tmp) > {log} 2>&1
        """


rule govaere_align:
    input: fq1=f"{GD}/fastq/{{run}}_1.fastq.gz", fq2=f"{GD}/fastq/{{run}}_2.fastq.gz",
           index=ancient(f"{REF}/star_v42_oh75/SA"), fa=ancient(GENOME), fai=ancient(GENOME + ".fai")
    output: cram=f"{GD}/cram/{{run}}.cram", crai=f"{GD}/cram/{{run}}.cram.crai",
            counts=f"{GD}/counts/{{run}}.ReadsPerGene.tsv", star_log=f"{GD}/star_logs/{{run}}.Log.final.out",
            sj=f"{GD}/star_logs/{{run}}.SJ.out.tab", ok=f"{GD}/cram/{{run}}.ok"
    threads: 8
    resources: mem_mb=64000, runtime=240     # 48 GB ran out of memory on some runs
    priority: 10          # align what is on disk before downloading more
    log: f"{GD}/logs/align/{{run}}.log"
    shell:
        "bash {SCRIPTS}/govaere_align.sh {wildcards.run} {input.fq1} {input.fq2} {REF}/star_v42_oh75 {input.fa} "
        "{GD} {threads} '{GOV[modules]}' > {log} 2>&1"
