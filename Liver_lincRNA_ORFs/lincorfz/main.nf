#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    LINCORFZ MAIN PIPELINE
    De novo annotation of expressed lincRNAs and prioritization of their translated ORFs
    from bulk RNA-seq, genomic Ribo-seq and polysome-fraction RNA-seq.
========================================================================================
*/

// main chromosomes (chrM is never analysed) or the user-defined subset
def chromList() {
    if (params.chromosomes && params.chromosomes.toString() != 'all') {
        return params.chromosomes.toString().tokenize(', ').collect { c -> c.trim() }.findAll { c -> c && c != 'chrM' }
    }
    def wanted = (1..22).collect { i -> "chr${i}".toString() } + ["chrX", "chrY"]
    def present = file(params.chrom_sizes).readLines().collect { line -> line.split('\t')[0] } as Set
    return wanted.findAll { c -> present.contains(c) }
}

def optionalFile(path, name) {
    return path ? file(path) : file("${projectDir}/assets/${name}")
}

workflow {

    ['rnaseq_tsv', 'riboseq_tsv', 'polysome_tsv', 'genome', 'chrom_sizes', 'reference_gtf', 'conservation_dir'].each { p ->
        if (!params[p]) { error("Please provide the parameter --${p}") }
    }

    def chroms     = chromList()
    def chrom_arg  = chroms.join(',')
    def restricted = params.chromosomes && params.chromosomes.toString() != 'all'
    def ribo_labels = file(params.riboseq_tsv).readLines().drop(1).findAll { line -> line.trim() }.collect { line -> line.split('\t')[0] }

    genome      = file(params.genome)
    genome_fai  = file("${params.genome}.fai")
    chrom_sizes = file(params.chrom_sizes)
    ref_gtf     = file(params.reference_gtf)
    rnaseq_tsv  = file(params.rnaseq_tsv)

    rnaseq_ch = channel.fromPath(params.rnaseq_tsv)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.sample, file(r.bam), file("${r.bam}.bai")) }
    ribo_ch = channel.fromPath(params.riboseq_tsv)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.sample, file(r.bam), file("${r.bam}.bai")) }

    // ==============================================================================
    // 1. Splice junctions (RNA-seq + annotation) -> whitelist for spliced Ribo-seq reads
    // ==============================================================================
    build_junctions(rnaseq_tsv, ref_gtf)

    // ==============================================================================
    // 2. Polysome profiling: per-nucleotide mean ribosome load (MRL) bigwigs
    // ==============================================================================
    mrl_in = channel.fromList(chroms).combine(channel.of('plus', 'minus'))
    polysome_mrl(mrl_in, file(params.polysome_tsv), chrom_sizes)
    polysome_bigwigs(
        polysome_mrl.out.bedgraphs.groupTuple().map { strand, fl -> tuple(strand, fl.flatten()) },
        chrom_sizes
    )

    // ==============================================================================
    // 3. Genomic Ribo-seq: P-site offset calibration and P-site tracks
    // ==============================================================================
    ribo_calibrate(ribo_ch, ref_gtf, build_junctions.out.whitelist, chrom_sizes, chrom_arg)
    ribo_psites(
        ribo_ch.join(ribo_calibrate.out.offsets),
        build_junctions.out.whitelist, chrom_sizes, chrom_arg
    )
    ribo_bigwigs(ribo_psites.out.bedgraphs.collect(), ribo_psites.out.stats.collect(), chrom_sizes)
    // read-length-specific priors on the frame distribution (canonical CDSs vs 3'UTRs)
    frame_priors(ribo_bigwigs.out.bigwigs, ref_gtf, chrom_sizes, chrom_arg)

    // ==============================================================================
    // 4. De novo transcriptome assembly and lincRNA selection
    // ==============================================================================
    if (restricted) {
        rnaseq_subset(rnaseq_ch, chrom_arg)
    }
    rna_bams = restricted ? rnaseq_subset.out.bam : rnaseq_ch
    stringtie_assemble(rna_bams, ref_gtf)
    stringtie_merge(stringtie_assemble.out.gtf.collect(), ref_gtf)
    gffcompare(stringtie_merge.out.gtf, ref_gtf)
    stringtie_quant(rna_bams, stringtie_merge.out.gtf)
    select_lincrnas(
        stringtie_merge.out.gtf, stringtie_quant.out.gtf.collect(), rnaseq_tsv, ref_gtf,
        build_junctions.out.junctions, optionalFile(params.noncode_gtf, 'NO_FILE_noncode'),
        chrom_sizes, chrom_arg
    )

    // ==============================================================================
    // 5. Expression of lincRNA loci (featureCounts + Sanity)
    // ==============================================================================
    featurecounts(
        rna_bams.map { _s, b, _i -> b }.collect(), rna_bams.map { _s, _b, i -> i }.collect(),
        ref_gtf, select_lincrnas.out.gtf
    )
    sanity(featurecounts.out.reference, featurecounts.out.lincrna, rnaseq_tsv)

    // ==============================================================================
    // 6. ORF calling and evidence layers
    // ==============================================================================
    call_orfs(select_lincrnas.out.gtf, select_lincrnas.out.tsv, genome, genome_fai, ribo_bigwigs.out.bigwigs)
    bayes_frame_test(call_orfs.out.length_counts, frame_priors.out.priors)
    polysome_scores(
        polysome_bigwigs.out.bigwigs.collect(), select_lincrnas.out.gtf, call_orfs.out.orfs,
        ref_gtf, chrom_sizes, chrom_arg
    )
    conservation(call_orfs.out.orfs, ref_gtf, chrom_sizes)

    if (params.peptide_tables) {
        peptide_ch = channel.fromList(params.peptide_tables.toString().tokenize(',').collect { f -> file(f.trim()) }).collect()
        canonical_proteome(genome, genome_fai, ref_gtf)
        match_peptides(call_orfs.out.orfs, peptide_ch, canonical_proteome.out.faa)
        peptides_out = match_peptides.out.tsv
    } else {
        peptides_out = channel.value(file("${projectDir}/assets/NO_FILE_peptides"))
    }

    // ==============================================================================
    // 7. Integration, ranking, browser files
    // ==============================================================================
    integrate_and_rank(
        call_orfs.out.orfs, select_lincrnas.out.tsv, select_lincrnas.out.gtf, call_orfs.out.tx_ribo,
        polysome_scores.out.transcripts, polysome_scores.out.orfs, sanity.out.lincrna,
        conservation.out.tsv, bayes_frame_test.out.tsv, peptides_out, ref_gtf,
        optionalFile(params.hcc_orfs_gtf, 'NO_FILE_hcc'), ribo_labels.join(',')
    )
    igv_session(integrate_and_rank.out.gtf, ribo_bigwigs.out.bigwigs, polysome_bigwigs.out.bigwigs.collect())
}

/*
========================================================================================
    PROCESS DEFINITIONS
========================================================================================
*/

process build_junctions {
    label 'process_low'
    label 'env_py'
    publishDir "${params.outdir}/junctions", mode: params.publish_mode

    input:
    path rnaseq_tsv
    path ref_gtf

    output:
    path "junctions.tsv.gz", emit: junctions
    path "whitelist.tsv.gz", emit: whitelist

    script:
    """
    build_junction_whitelist.py \\
        --rnaseq-tsv ${rnaseq_tsv} \\
        --gtf ${ref_gtf} \\
        --min-unique-reads ${params.junction_min_unique_reads}
    """
}

process polysome_mrl {
    tag "${chrom} ${strand}"
    label 'process_low'
    label 'env_py'

    input:
    tuple val(chrom), val(strand)
    path polysome_tsv
    path chrom_sizes

    output:
    tuple val(strand), path("${chrom}.${strand}.*.bedGraph"), emit: bedgraphs

    script:
    def per_rep = params.mrl_per_replicate ? '--per-replicate' : ''
    """
    compute_mrl.py \\
        --table ${polysome_tsv} \\
        --chrom ${chrom} \\
        --strand ${strand} \\
        --chrom-sizes ${chrom_sizes} \\
        --min-ribosomes ${params.mrl_min_ribosomes} \\
        --max-ribosomes ${params.mrl_max_ribosomes} \\
        --min-total-coverage ${params.mrl_min_total_coverage} \\
        ${per_rep}
    """
}

process polysome_bigwigs {
    tag "${strand}"
    label 'process_medium_low_cpu'
    label 'env_tools'
    publishDir "${params.outdir}/polysome_mrl", mode: params.publish_mode

    input:
    tuple val(strand), path(bedgraphs)
    path chrom_sizes

    output:
    path "polysome.*.${strand}.bw", emit: bigwigs

    script:
    """
    for track in \$(ls *.${strand}.*.bedGraph | sed -E 's/^[^.]+\\.${strand}\\.//; s/\\.bedGraph\$//' | sort -u); do
        cat *.${strand}.\${track}.bedGraph | LC_ALL=C sort -k1,1 -k2,2n -S 4G > merged.bg
        if [ -s merged.bg ]; then
            bedGraphToBigWig merged.bg ${chrom_sizes} polysome.\${track}.${strand}.bw
        else
            echo "[warning] no data for track \${track} (${strand})" >&2
        fi
        rm merged.bg
    done
    """
}

process ribo_calibrate {
    tag "${sample}"
    label 'process_medium_low_cpu'
    label 'env_py'
    publishDir "${params.outdir}/riboseq/calibration", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_gtf
    path whitelist
    path chrom_sizes
    val chrom_arg

    output:
    tuple val(sample), path("${sample}.offsets.tsv"), emit: offsets
    path "${sample}.metagene.tsv", emit: metagene
    path "${sample}.calibration.pdf", emit: pdf

    script:
    """
    riboseq_calibrate_offsets.py \\
        --bam ${bam} \\
        --gtf ${ref_gtf} \\
        --whitelist ${whitelist} \\
        --chrom-sizes ${chrom_sizes} \\
        --chromosomes ${chrom_arg} \\
        --prefix ${sample} \\
        --min-len ${params.ribo_min_len} \\
        --max-len ${params.ribo_max_len} \\
        --max-clip5 ${params.ribo_max_clip5} \\
        --min-offset ${params.ribo_min_offset} \\
        --max-offset ${params.ribo_max_offset} \\
        --min-frame0 ${params.ribo_min_frame0}
    """
}

process ribo_psites {
    tag "${sample}"
    label 'process_medium_low_cpu'
    label 'env_py'
    publishDir "${params.outdir}/riboseq/psite_stats", mode: params.publish_mode, pattern: "*.psite_stats.tsv"

    input:
    tuple val(sample), path(bam), path(bai), path(offsets)
    path whitelist
    path chrom_sizes
    val chrom_arg

    output:
    path "${sample}.*.bedGraph", emit: bedgraphs
    path "${sample}.psite_stats.tsv", emit: stats

    script:
    """
    riboseq_psites.py \\
        --bam ${bam} \\
        --offsets ${offsets} \\
        --whitelist ${whitelist} \\
        --chrom-sizes ${chrom_sizes} \\
        --chromosomes ${chrom_arg} \\
        --max-clip5 ${params.ribo_max_clip5} \\
        --min-len ${params.ribo_min_len} \\
        --max-len ${params.ribo_max_len} \\
        --min-cds-reads-length-tracks ${params.ribo_length_min_cds_reads} \\
        --prefix ${sample}
    """
}

process ribo_bigwigs {
    label 'process_medium_low_cpu'
    label 'env_tools'
    // main tracks -> psite_bigwigs/; per-read-length CPM tracks -> psite_bigwigs/by_read_length/;
    // per-read-length raw counts are only used downstream (frame priors, ORF frame counts)
    publishDir "${params.outdir}/riboseq/psite_bigwigs", mode: params.publish_mode,
        saveAs: { fn -> fn.endsWith('.raw.bw') ? null : (fn.contains('.len') ? "by_read_length/${fn}" : fn) }

    input:
    path bedgraphs
    path stats
    path chrom_sizes

    output:
    path "*.bw", emit: bigwigs

    script:
    """
    samples=\$(ls *.psite_stats.tsv | sed 's/\\.psite_stats\\.tsv\$//')
    total=0
    total_len=0
    for s in \$samples; do
        n=\$(awk -F'\\t' '\$1=="psites"{print \$2}' \${s}.psite_stats.tsv)
        m=\$(awk -F'\\t' '\$1=="psites_all_lengths"{print \$2}' \${s}.psite_stats.tsv)
        total=\$((total + n))
        total_len=\$((total_len + m))
    done
    echo "pooled P-sites (main tracks): \$total; all read lengths: \$total_len"

    to_cpm() {  # bedGraph, total, output bigwig
        awk -F'\\t' -v t=\$2 'BEGIN{OFS="\\t"} {printf "%s\\t%s\\t%s\\t%.6f\\n", \$1,\$2,\$3,\$4*1e6/t}' \$1 > \$1.cpm.tmp
        bedGraphToBigWig \$1.cpm.tmp ${chrom_sizes} \$3
        rm \$1.cpm.tmp
    }
    sum_union() {  # output bedGraph, inputs...
        out=\$1; shift
        if [ \$# -gt 1 ]; then
            bedtools unionbedg -i "\$@" \\
                | awk -F'\\t' 'BEGIN{OFS="\\t"} {s=0; for(i=4;i<=NF;i++) s+=\$i; print \$1,\$2,\$3,s}' > \$out
        else
            cp \$1 \$out
        fi
    }

    # ---- main tracks (read classes that passed the calibration)
    for strand in plus minus; do
        for kind in "" "unique."; do
            sorted=()
            for s in \$samples; do
                LC_ALL=C sort -k1,1 -k2,2n -S 2G \${s}.\${kind}\${strand}.bedGraph > \${s}.\${kind}\${strand}.sorted.bg
                sorted+=(\${s}.\${kind}\${strand}.sorted.bg)
                if [ -z "\${kind}" ] && [ -s \${s}.\${strand}.sorted.bg ]; then
                    bedGraphToBigWig \${s}.\${strand}.sorted.bg ${chrom_sizes} \${s}.\${strand}.bw
                fi
            done
            label=\$([ -z "\${kind}" ] && echo pooled || echo pooled_unique)
            sum_union \${label}.\${strand}.bg \${sorted[@]}
            bedGraphToBigWig \${label}.\${strand}.bg ${chrom_sizes} \${label}.\${strand}.bw
            if [ "\${label}" = "pooled" ]; then
                to_cpm pooled.\${strand}.bg \$total pooled_cpm.\${strand}.bw
            fi
            rm -f \${sorted[@]}
        done
    done

    # ---- per aligned read length: raw counts and CPM (scaled by the sample's P-sites over all lengths)
    for L in \$(ls *.len*.plus.bedGraph | sed -E 's/.*\\.len([0-9]+)\\.plus\\.bedGraph/\\1/' | sort -un); do
        for strand in plus minus; do
            sorted=()
            for s in \$samples; do
                f=\${s}.len\${L}.\${strand}.bedGraph
                [ -s \$f ] || continue
                tot=\$(awk -F'\\t' '\$1=="psites_all_lengths"{print \$2}' \${s}.psite_stats.tsv)
                LC_ALL=C sort -k1,1 -k2,2n -S 2G \$f > \${s}.len\${L}.\${strand}.sorted.bg
                bedGraphToBigWig \${s}.len\${L}.\${strand}.sorted.bg ${chrom_sizes} \${s}.len\${L}.\${strand}.raw.bw
                to_cpm \${s}.len\${L}.\${strand}.sorted.bg \$tot \${s}.len\${L}.\${strand}.cpm.bw
                sorted+=(\${s}.len\${L}.\${strand}.sorted.bg)
            done
            if [ \${#sorted[@]} -gt 0 ]; then
                sum_union pooled.len\${L}.\${strand}.bg \${sorted[@]}
                to_cpm pooled.len\${L}.\${strand}.bg \$total_len pooled.len\${L}.\${strand}.cpm.bw
            fi
            rm -f \${sorted[@]}
        done
    done
    rm -f *.bg
    """
}

process frame_priors {
    label 'process_medium'
    label 'env_py'
    publishDir "${params.outdir}/riboseq/frame_priors", mode: params.publish_mode

    input:
    path bigwigs
    path ref_gtf
    path chrom_sizes
    val chrom_arg

    output:
    path "frame_priors.tsv", emit: priors
    path "frame_counts_by_length.tsv"
    path "frame_priors.pdf"

    script:
    """
    printf "sample\\tlength\\tstrand\\tpath\\n" > length_tracks.tsv
    for f in *.len*.plus.raw.bw; do
        s=\${f%%.len*}
        L=\$(echo \$f | sed -E 's/.*\\.len([0-9]+)\\.plus\\.raw\\.bw/\\1/')
        printf "%s\\t%s\\t+\\t%s\\n" "\$s" "\$L" "\$f" >> length_tracks.tsv
        [ -e \${s}.len\${L}.minus.raw.bw ] && printf "%s\\t%s\\t-\\t%s\\n" "\$s" "\$L" "\${s}.len\${L}.minus.raw.bw" >> length_tracks.tsv
    done
    frame_priors.py \\
        --length-tracks length_tracks.tsv \\
        --gtf ${ref_gtf} \\
        --chrom-sizes ${chrom_sizes} \\
        --chromosomes ${chrom_arg} \\
        --position-cap ${params.bayes_count_positions ? 1 : 0} \\
        --window-codons ${params.bayes_window_codons} \\
        --n-cores ${task.cpus}
    """
}

process bayes_frame_test {
    label 'process_medium_low_cpu'
    label 'env_py'
    publishDir "${params.outdir}/orfs", mode: params.publish_mode

    input:
    path counts
    path priors

    output:
    path "orf_frame_bayes.tsv", emit: tsv

    script:
    """
    bayes_frame_test.py \\
        --counts ${counts} \\
        --priors ${priors} \\
        --prior-h1 ${params.bayes_prior_h1} \\
        --w-noise ${params.bayes_w_noise} \\
        --use ${params.bayes_count_positions ? 'positions' : 'reads'}
    """
}

process rnaseq_subset {
    tag "${sample}"
    label 'process_low'
    label 'env_tools'

    input:
    tuple val(sample), path(bam), path(bai)
    val chrom_arg

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam

    script:
    """
    samtools view -b -o ${sample}.bam ${bam} ${chrom_arg.tokenize(',').join(' ')}
    samtools index ${sample}.bam
    """
}

process stringtie_assemble {
    tag "${sample}"
    label 'process_medium'
    label 'env_tools'
    publishDir "${params.outdir}/assembly/per_sample", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_gtf

    output:
    path "${sample}.stringtie.gtf", emit: gtf

    script:
    """
    stringtie ${bam} --rf -G ${ref_gtf} -p ${task.cpus} -l ${sample} ${params.stringtie_args} -o ${sample}.stringtie.gtf
    """
}

process stringtie_merge {
    label 'process_medium'
    label 'env_tools'
    publishDir "${params.outdir}/assembly", mode: params.publish_mode

    input:
    path gtfs
    path ref_gtf

    output:
    path "merged.gtf", emit: gtf

    script:
    """
    stringtie --merge -G ${ref_gtf} ${params.stringtie_merge_args} -o merged.gtf ${gtfs}
    """
}

process gffcompare {
    label 'process_low'
    label 'env_tools'
    publishDir "${params.outdir}/assembly/gffcompare", mode: params.publish_mode

    input:
    path merged
    path ref_gtf

    output:
    path "gffcmp.*"

    script:
    """
    gffcompare -r ${ref_gtf} -o gffcmp ${merged}
    """
}

process stringtie_quant {
    tag "${sample}"
    label 'process_medium'
    label 'env_tools'
    publishDir "${params.outdir}/assembly/quantification", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    path merged

    output:
    path "${sample}.quant.gtf", emit: gtf
    path "${sample}.gene_abund.tab", emit: genes

    script:
    """
    stringtie ${bam} --rf -e -G ${merged} -p ${task.cpus} -A ${sample}.gene_abund.tab -o ${sample}.quant.gtf
    """
}

process select_lincrnas {
    label 'process_medium_low_cpu'
    label 'env_py'
    publishDir "${params.outdir}/lincRNAs", mode: params.publish_mode

    input:
    path merged
    path quant_gtfs
    path rnaseq_tsv
    path ref_gtf
    path junctions
    path noncode_gtf
    path chrom_sizes
    val chrom_arg

    output:
    path "lincRNA_transcripts.gtf", emit: gtf
    path "lincRNA_transcripts.tsv", emit: tsv
    path "lincRNA_candidates.all.tsv", emit: all

    script:
    def nc = noncode_gtf.name.startsWith('NO_FILE') ? '' : "--noncode-gtf ${noncode_gtf}"
    def pc_exon = params.keep_pc_gene_exon_overlap ? '--keep-pc-gene-exon-overlap' : ''
    def pc_intr = params.keep_pc_intronic ? '--keep-pc-intronic' : ''
    """
    select_lincrnas.py \\
        --merged-gtf ${merged} \\
        --quant-gtfs ${quant_gtfs} \\
        --rnaseq-tsv ${rnaseq_tsv} \\
        --reference-gtf ${ref_gtf} \\
        --junctions ${junctions} \\
        ${nc} \\
        --chrom-sizes ${chrom_sizes} \\
        --chromosomes ${chrom_arg} \\
        --min-tx-length ${params.min_tx_length} \\
        --min-tpm ${params.min_tpm} \\
        --min-tpm-monoexonic ${params.min_tpm_monoexonic} \\
        --min-junction-reads ${params.min_junction_reads} \\
        --min-isoform-fraction ${params.min_isoform_fraction} \\
        --utr-extension-distance ${params.utr_extension_distance} \\
        ${pc_exon} ${pc_intr}
    """
}

process featurecounts {
    label 'process_medium'
    label 'env_tools'
    publishDir "${params.outdir}/expression", mode: params.publish_mode

    input:
    path bams
    path bais
    path ref_gtf
    path lincrna_gtf

    output:
    path "counts.reference.txt", emit: reference
    path "counts.lincRNA_loci.txt", emit: lincrna
    path "*.summary"

    script:
    """
    featureCounts -p --countReadPairs -s ${params.featurecounts_strand} -T ${task.cpus} \\
        -a ${ref_gtf} -t exon -g gene_id -o counts.reference.txt ${bams}
    featureCounts -p --countReadPairs -s ${params.featurecounts_strand} -T ${task.cpus} -O \\
        -a ${lincrna_gtf} -t exon -g gene_id -o counts.lincRNA_loci.txt ${bams}
    """
}

process sanity {
    label 'process_high'
    label 'env_py'
    publishDir "${params.outdir}/expression", mode: params.publish_mode

    input:
    path ref_counts
    path linc_counts
    path rnaseq_tsv

    output:
    path "lincRNA_loci.expression.tsv", emit: lincrna
    path "sanity.*.tsv"

    script:
    """
    export PYTHONNOUSERSITE=1
    run_sanity.py \\
        --reference-counts ${ref_counts} \\
        --lincrna-counts ${linc_counts} \\
        --rnaseq-tsv ${rnaseq_tsv} \\
        --min-total-counts ${params.sanity_min_total_counts} \\
        --v-method ${params.sanity_v_method} \\
        --n-threads ${task.cpus}
    """
}

process call_orfs {
    label 'process_high'
    label 'env_py'
    publishDir "${params.outdir}/orfs", mode: params.publish_mode

    input:
    path lincrna_gtf
    path lincrna_tsv
    path genome
    path genome_fai
    path bigwigs

    output:
    path "orfs.tsv", emit: orfs
    path "orf_isoforms.tsv", emit: isoforms
    path "transcript_ribo.tsv", emit: tx_ribo
    path "orf_length_frame_counts.tsv", emit: length_counts
    path "psite_tracks.tsv"
    path "length_tracks.tsv"

    script:
    """
    printf "label\\tstrand\\tpath\\n" > psite_tracks.tsv
    for f in *.plus.bw; do
        lab=\${f%.plus.bw}
        [ "\$lab" = "pooled_cpm" ] && continue
        printf "%s\\t+\\t%s\\n" "\$lab" "\$f" >> psite_tracks.tsv
        printf "%s\\t-\\t%s\\n" "\$lab" "\${lab}.minus.bw" >> psite_tracks.tsv
    done
    printf "sample\\tlength\\tstrand\\tpath\\n" > length_tracks.tsv
    for f in *.len*.plus.raw.bw; do
        s=\${f%%.len*}
        L=\$(echo \$f | sed -E 's/.*\\.len([0-9]+)\\.plus\\.raw\\.bw/\\1/')
        printf "%s\\t%s\\t+\\t%s\\n" "\$s" "\$L" "\$f" >> length_tracks.tsv
        [ -e \${s}.len\${L}.minus.raw.bw ] && printf "%s\\t%s\\t-\\t%s\\n" "\$s" "\$L" "\${s}.len\${L}.minus.raw.bw" >> length_tracks.tsv
    done
    call_orfs.py \\
        --length-tracks length_tracks.tsv \\
        --transcripts-gtf ${lincrna_gtf} \\
        --transcripts-tsv ${lincrna_tsv} \\
        --genome ${genome} \\
        --psite-tracks psite_tracks.tsv \\
        --min-aa ${params.min_aa} \\
        --start-codons ${params.start_codons} \\
        --near-cognate-codons "${params.near_cognate_codons}" \\
        --min-psites ${params.min_psites_orf} \\
        --n-cores ${task.cpus}
    """
}

process polysome_scores {
    label 'process_medium'
    label 'env_py'
    publishDir "${params.outdir}/polysome_scores", mode: params.publish_mode

    input:
    path bigwigs
    path lincrna_gtf
    path orfs
    path ref_gtf
    path chrom_sizes
    val chrom_arg

    output:
    path "polysome.transcripts.tsv", emit: transcripts
    path "polysome.orfs.tsv", emit: orfs
    path "polysome.reference_mRNAs.tsv"
    path "polysome_mrl.pdf"

    script:
    """
    printf "n_ribosomes\\tstrand\\tpath\\n" > fraction_tracks.tsv
    for f in polysome.frac*.plus.bw; do
        k=\$(echo \$f | sed -E 's/polysome\\.frac([0-9]+)\\.plus\\.bw/\\1/')
        printf "%s\\t+\\t%s\\n" "\$k" "\$f" >> fraction_tracks.tsv
        printf "%s\\t-\\t%s\\n" "\$k" "polysome.frac\${k}.minus.bw" >> fraction_tracks.tsv
    done
    polysome_scores.py \\
        --fraction-tracks fraction_tracks.tsv \\
        --transcripts-gtf ${lincrna_gtf} \\
        --orfs ${orfs} \\
        --reference-gtf ${ref_gtf} \\
        --chrom-sizes ${chrom_sizes} \\
        --chromosomes ${chrom_arg} \\
        --min-mean-coverage ${params.polysome_min_mean_coverage} \\
        --n-cores ${task.cpus}
    """
}

process conservation {
    label 'process_medium_low_cpu'
    label 'env_py'
    publishDir "${params.outdir}/conservation", mode: params.publish_mode

    input:
    path orfs
    path ref_gtf
    path chrom_sizes

    output:
    path "conservation.orfs.tsv", emit: tsv

    script:
    """
    conservation_scores.py \\
        --orfs ${orfs} \\
        --reference-gtf ${ref_gtf} \\
        --conservation-dir ${params.conservation_dir} \\
        --chrom-sizes ${chrom_sizes}
    """
}

process canonical_proteome {
    label 'process_low'
    label 'env_tools'

    input:
    path genome
    path genome_fai
    path ref_gtf

    output:
    path "canonical_proteome.faa", emit: faa

    script:
    """
    gffread ${ref_gtf} -g ${genome} -y canonical_proteome.faa
    """
}

process match_peptides {
    label 'process_medium_low_cpu'
    label 'env_py'
    publishDir "${params.outdir}/peptides", mode: params.publish_mode

    input:
    path orfs
    path tables
    path proteome

    output:
    path "peptides.orfs.tsv", emit: tsv

    script:
    """
    match_peptides.py --orfs ${orfs} --peptide-tables ${tables} --proteome-fasta ${proteome}
    """
}

process integrate_and_rank {
    label 'process_medium_low_cpu'
    label 'env_py'
    publishDir "${params.outdir}/final", mode: params.publish_mode

    input:
    path orfs
    path lincrna_tsv
    path lincrna_gtf
    path tx_ribo
    path poly_tx
    path poly_orfs
    path expression
    path cons
    path bayes
    path peptides
    path ref_gtf
    path hcc_gtf
    val replicate_labels

    output:
    path "lincORFs.igv.gtf", emit: gtf
    path "lincORFs.*"

    script:
    def pep = peptides.name.startsWith('NO_FILE') ? '' : "--peptides ${peptides}"
    def hcc = hcc_gtf.name.startsWith('NO_FILE') ? '' : "--hcc-orfs-gtf ${hcc_gtf}"
    """
    integrate_and_rank.py \\
        --orfs ${orfs} \\
        --transcripts-tsv ${lincrna_tsv} \\
        --transcripts-gtf ${lincrna_gtf} \\
        --transcript-ribo ${tx_ribo} \\
        --polysome-transcripts ${poly_tx} \\
        --polysome-orfs ${poly_orfs} \\
        --expression ${expression} \\
        --conservation ${cons} \\
        --bayes ${bayes} \\
        --bayes-strong ${params.bayes_strong} \\
        --bayes-moderate ${params.bayes_moderate} \\
        --bayes-weak ${params.bayes_weak} \\
        --bayes-rep-support ${params.bayes_rep_support} \\
        ${pep} ${hcc} \\
        --reference-gtf ${ref_gtf} \\
        --replicate-labels ${replicate_labels} \\
        --primary-condition ${params.primary_condition} \\
        --p-strong ${params.p_strong} \\
        --min-inframe-strong ${params.min_inframe_strong} \\
        --min-frac-inframe-strong ${params.min_frac_inframe_strong} \\
        --min-reps-strong ${params.min_reps_strong} \\
        --p-moderate ${params.p_moderate} \\
        --min-inframe-moderate ${params.min_inframe_moderate} \\
        --p-weak ${params.p_weak} \\
        --min-frac-unique ${params.min_frac_unique} \\
        --frame-test-min-codons ${params.frame_test_min_codons} \\
        --frame-test-min-frac-inframe ${params.frame_test_min_frac_inframe} \\
        ${params.frame_test ? '--frame-test' : ''} \\
        --w-mrl ${params.w_mrl} --w-ribo ${params.w_ribo} --w-island ${params.w_island} \\
        --w-expr ${params.w_expr} --w-cons ${params.w_cons} --w-unique ${params.w_unique} \\
        --w-pep ${params.w_pep} \\
        --gtf-min-evidence ${params.gtf_min_evidence}
    """
}

process igv_session {
    label 'process_low'
    label 'env_py'
    publishDir "${params.outdir}/final", mode: params.publish_mode

    input:
    path final_gtf
    path psite_bws
    path polysome_bws

    output:
    path "lincorfz.igv_session.xml"
    path "igv_tracks.tsv"

    script:
    def out  = file(params.outdir).toAbsolutePath().toString()
    def cons = file(params.conservation_dir).toAbsolutePath().toString()
    def rows = []
    rows << ['feature', '-', 'lincorfz: lincRNAs + ORFs', "${out}/final/lincORFs.igv.gtf", '0,0,0', '-', '60']
    rows << ['feature', '-', 'lincorfz: ORFs (BED, colour = evidence)', "${out}/final/lincORFs.orfs.bed", '0,0,0', '-', '40']
    if (params.igv_reference_annotation) {
        rows << ['feature', '-', 'reference annotation', params.igv_reference_annotation, '0,0,0', '-', '60']
    }
    rows << ['data', 'Mean ribosome load (80S..decasome)', 'MRL +', "${out}/polysome_mrl/polysome.mrl.plus.bw", '0,102,204', '1,10', '60']
    rows << ['data', 'Mean ribosome load (80S..decasome)', 'MRL -', "${out}/polysome_mrl/polysome.mrl.minus.bw", '204,0,0', '1,10', '60']
    rows << ['data', 'Polysomal coverage (sum of fractions)', 'polysomal coverage +', "${out}/polysome_mrl/polysome.total.plus.bw", '0,102,204', '-', '50']
    rows << ['data', 'Polysomal coverage (sum of fractions)', 'polysomal coverage -', "${out}/polysome_mrl/polysome.total.minus.bw", '204,0,0', '-', '50']
    rows << ['data', 'Ribo-seq P-sites (CPM, pooled)', 'P-sites +', "${out}/riboseq/psite_bigwigs/pooled_cpm.plus.bw", '0,102,204', '-', '60']
    rows << ['data', 'Ribo-seq P-sites (CPM, pooled)', 'P-sites -', "${out}/riboseq/psite_bigwigs/pooled_cpm.minus.bw", '204,0,0', '-', '60']
    if (params.rnaseq_tracks_tsv) {
        file(params.rnaseq_tracks_tsv).readLines().drop(1).findAll { line -> line.trim() }.each { line ->
            def f = line.split('\t')
            rows << ['data', "RNA-seq ${f[0]}", "${f[0]} ${f[1]}", f[2], f[1] == '+' ? '0,102,204' : '204,0,0', '-', '50']
        }
    }
    [1, 2, 3].each { fr -> rows << ['data', 'PhyloCSF + strand', "PhyloCSF+${fr}", "${cons}/PhyloCSF+${fr}.bw", ['255,51,0', '0,0,204', '102,204,0'][fr - 1], '-', '40'] }
    [1, 2, 3].each { fr -> rows << ['data', 'PhyloCSF - strand', "PhyloCSF-${fr}", "${cons}/PhyloCSF-${fr}.bw", ['255,51,204', '0,153,204', '153,153,0'][fr - 1], '-', '40'] }
    rows << ['data', '-', 'PhyloCSF_Power', "${cons}/PhyloCSF_Power.bw", '100,100,100', '-', '20']
    def table = (['panel\tgroup\tname\tpath\tcolor\trange\theight'] + rows.collect { r -> r.join('\t') }).join('\n')
    """
    cat > igv_tracks.tsv <<'EOF'
${table}
EOF
    # P-sites per aligned read length (pooled samples, CPM), one colour per length
    palette=(31,119,180 255,127,14 44,160,44 214,39,40 148,103,189 140,86,75 227,119,194 127,127,127 188,189,34 23,190,207 0,0,128 128,0,0 0,128,128 128,128,0)
    for strand in plus minus; do
        sign=\$([ \$strand = plus ] && echo "+" || echo "-")
        i=0
        for L in \$(ls pooled.len*.\${strand}.cpm.bw 2>/dev/null | sed -E 's/pooled\\.len([0-9]+)\\..*/\\1/' | sort -n); do
            printf "data\\tRibo-seq P-sites by read length (%s strand, pooled CPM)\\tlength %s nt %s\\t%s\\t%s\\t-\\t70\\n" \\
                "\$sign" "\$L" "\$sign" "${out}/riboseq/psite_bigwigs/by_read_length/pooled.len\${L}.\${strand}.cpm.bw" \\
                "\${palette[\$((i % \${#palette[@]}))]}" >> igv_tracks.tsv
            i=\$((i + 1))
        done
    done
    make_igv_session.py --genome ${file(params.genome).toAbsolutePath()} --tracks igv_tracks.tsv \\
        --locus ${params.igv_locus} --out lincorfz.igv_session.xml
    """
}
