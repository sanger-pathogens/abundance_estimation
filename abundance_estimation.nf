#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { metawrap_qc } from './modules/metawrap_qc.nf'
include { cleanup_sorted_bam_files; cleanup_trimmed_fastq_files; cleanup_instrain_output } from './modules/cleanup.nf'

// helper functions
def printHelp() {
    log.info """
    Usage:
    nextflow run abundance_estimation.nf

    Options:
      --manifest                   Manifest containing paths to fastq files (mandatory)
      --btidx                      bowtie index - default: /lustre/scratch118/infgen/team162/shared/gtdb_genomes_reps_r202/gtdb_genomes_reps_r202.fasta.bt2 (optional)
      --genome_file                genome file - default: /lustre/scratch118/infgen/team162/shared/gtdb_genomes_reps_r202/gtdb_genomes_reps_r202.fasta (optional)
      --stb_file                   stb file - default: /lustre/scratch118/infgen/team162/shared/gtdb_genomes_reps_r202/gtdb_genomes_reps_r202.stb
      --bowtie2_samtools_threads   threads - default: 16 (optional)
      --instrain_threads           threads - default: 16 (optional)
      --instrain_full_output       get full instrain output - default false (optional)
      --keep_metawrap_qc           don't cleanup metawrap_qc output - default false (optional)
      --keep_bowtie2samtools       don't cleanup bowtie2samtools output - default false (optional)
      --keep_instrain              don't cleanup instrain output - default false (optional)
      --skip_qc                    skip metawrap qc step - default false (optional)
      --instrain_queue             job queue for instrain - default normal (optional)
      --bowtie2samtools_queue      job queue for bowtie2samtools - default normal (optional)
      -profile                     always use sanger_lsf when running on the farm (mandatory)
      --help                       print this help message (optional)
    """.stripIndent()
}


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
    path("${genome_info_file}"), optional: true
    path(sorted_bam), emit: sorted_bam
    path("${workdir}"), emit: workdir
    val sample_id, emit: sample_id

    script:
    genome_info_file="${sample_id}_genome_info.tsv"
    workdir="workdir.txt"
    """
    pwd > workdir.txt
    inStrain profile $sorted_bam $genome_file -o $sample_id -p $threads -s $stb_file --database_mode --skip_plot_generation
    if $params.instrain_full_output
    then
        mkdir -p ${workflow.projectDir}/${params.results_dir}/${sample_id}
        cp -r $sample_id/* ${workflow.projectDir}/${params.results_dir}/${sample_id}
    else
        mv ${sample_id}/output/${sample_id}"_genome_info.tsv" .
    fi
    """
}

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }
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
    if (!params.keep_metawrap_qc) {
        if (!params.skip_qc) {
            cleanup_trimmed_fastq_files(bowtie2samtools.out.trimmed_fastqs)
        }
    }
    instrain(bowtie2samtools.out.bam_file, params.genome_file, params.stb_file, params.instrain_threads)
    if (!params.keep_bowtie2samtools) {
        cleanup_sorted_bam_files(instrain.out.sorted_bam)
    }
    if (!params.keep_instrain) {
        cleanup_instrain_output(instrain.out.workdir, instrain.out.sample_id)
    }
}
