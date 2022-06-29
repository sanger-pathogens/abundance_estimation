#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { metawrap_qc } from './modules/metawrap_qc.nf'
include { cleanup_sorted_bam_files } from './modules/cleanup.nf'
include { cleanup_trimmed_fastq_files } from './modules/cleanup.nf'
include { cleanup_instrain_output } from './modules/cleanup.nf'

def validate_parameters() {
    // Parameter checking function
    def errors = 0

    if (params.manifest) {
        manifest=file(params.manifest)
        if (!manifest.exists()) {
            log.error("The manifest file specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No manifest file specified. Please specify one using the --manifest option.")
        errors += 1
    }

    if (params.genome_file) {
        genome_file=file(params.genome_file)
        if (!genome_file.exists()) {
            log.error("The genome file specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No genome file specified. Please specify one using the --genome_file option.")
        errors += 1
    }

    if (params.stb_file) {
        stb_file=file(params.stb_file)
        if (!stb_file.exists()) {
            log.error("The stb file specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No stb file specified. Please specify one using the --stb_file option.")
        errors += 1
    }

    if (params.bowtie2_samtools_threads) {
        if (!params.bowtie2_samtools_threads.toString().isInteger()) {
        log.error("Please ensure the bowtie2_samtools_threads parameter is a number")
        errors += 1
        }
    }
    else {
        log.error("Please specify the number of threads for bowtie2samtools using the --bowtie2_samtools_threads option")
        errors += 1
    }

    if (params.instrain_threads) {
        if (!params.instrain_threads.toString().isInteger()) {
        log.error("Please ensure the instrain_threads parameter is a number")
        errors += 1
        }
    }
    else {
        log.error("Please specify the number of threads for instrain_threads using the --instrain_threads option")
        errors += 1
    }

    if (params.results_dir) {
        results_dir_path=file(params.results_dir)
        if (!results_dir_path.getParent().exists()) {
            log.error("The results directory path specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No results directory has been specified, please ensure you provide a value or a default.")
        errors += 1
    }

    if (errors > 0) {
            log.error(String.format("%d errors detected", errors))
            exit 1
        }
}


process bowtie2samtools {
    input:
    tuple val(sample_id), file(first_read), file(second_read)
    val btidx
    val threads

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam_file
    tuple val(sample_id), path(first_read), path(second_read), emit: trimmed_fastqs

    script:
    """
    bowtie2 -p $threads -x $btidx -1 $first_read -2 $second_read | samtools sort -@ $threads -o ${sample_id}".sorted.bam"
    """
}

process instrain {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*.tsv"
    input:
    tuple val(sample_id), file(sorted_bam)
    val genome_file
    val stb_file
    val threads

    output:
    path("${genome_info_file}")
    path(sorted_bam), emit: sorted_bam
    path("${workdir}"), emit: workdir
    val sample_id, emit: sample_id

    script:
    genome_info_file="${sample_id}_genome_info.tsv"
    workdir="workdir.txt"
    """
    pwd > workdir.txt
    inStrain profile $sorted_bam $genome_file -o $sample_id -p $threads -s $stb_file --database_mode --skip_plot_generation
    mv ${sample_id}/output/${sample_id}"_genome_info.tsv" .
    """
}

process instrain_full_output {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true
    input:
    tuple val(sample_id), file(sorted_bam)
    val genome_file
    val stb_file
    val threads

    output:
    path("${sample_id}/*")

    script:
    """
    inStrain profile $sorted_bam $genome_file -o $sample_id -p $threads -s $stb_file --database_mode --skip_plot_generation
    """
}

workflow {
    validate_parameters()
    manifest_ch = Channel.fromPath(params.manifest)
    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.first_read), file(row.second_read)) }
    if (params.skip_qc) {
        bowtie2samtools(fastq_path_ch, params.btidx, params.bowtie2_samtools_threads)
    }
    else {
        metawrap_qc(fastq_path_ch)
        bowtie2samtools(metawrap_qc.out.trimmed_fastqs, params.btidx, params.bowtie2_samtools_threads)
    }
    if (params.cleanup_metawrapqc) {
        if (!params.skip_qc) {
            cleanup_trimmed_fastq_files(bowtie2samtools.out.trimmed_fastqs)
        }
    }
    if (params.full_output) {
        // presuming if you are using full output, it is for debugging purposes, so no clean up
        instrain_full_output(bowtie2samtools.out.bam_file, params.genome_file, params.stb_file, params.instrain_threads)
    }
    else {
        instrain(bowtie2samtools.out.bam_file, params.genome_file, params.stb_file, params.instrain_threads)
        if (params.cleanup_bowtie2samtools) {
            cleanup_sorted_bam_files(instrain.out.sorted_bam)
        }
        if (params.cleanup_instrain_output) {
            cleanup_instrain_output(instrain.out.workdir, instrain.out.sample_id)
        }
    }
}
