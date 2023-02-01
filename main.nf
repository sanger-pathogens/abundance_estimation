#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { metawrap_qc } from './modules/metawrap_qc.nf'
include { cleanup_sorted_bam_files; cleanup_trimmed_fastq_files; cleanup_instrain_output } from './modules/cleanup.nf'

// helper functions
def printHelp() {
    log.info """
    Usage:
    nextflow run .

    Options:
      --manifest                   Manifest containing paths to fastq files (mandatory)
      --btidx                      bowtie index - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.bt2 (optional)
      --genome_file                genome file - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.fasta (optional)
      --stb_file                   stb file - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.stb (optional)
      --bowtie2_samtools_threads   threads - default: 8 (optional)
      --instrain_threads           threads - default: 8 (optional)
      --instrain_full_output       get full instrain output - default false (optional)
      --keep_metawrap_qc           don't cleanup metawrap_qc output - default false (optional)
      --keep_bowtie2samtools       don't cleanup bowtie2samtools output - default false (optional)
      --keep_instrain              don't cleanup instrain output - default false (optional)
      --skip_qc                    skip metawrap qc step - default false (optional)
      --instrain_queue             job queue for instrain - default long (optional)
      --instrain_quick_profile     use quick-profile option for inStrain - default false (optional)
      --bowtie2_samtools_only      only run bowtie2_samtools process - default false (optional)
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

    if (params.instrain_full_output && params.instrain_quick_profile) {
        log.error("the --instrain_full_output and --instrain_quick_profile options are incompatible, please choose one of these options.")
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
    container '/software/pathogen/images/bowtie2-samtools-1.0.simg'
    if (params.bowtie2_samtools_only) { publishDir path: "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*.sorted.bam" }
    input:
    tuple val(sample_id), file(first_read), file(second_read)
    val btidx
    val threads

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam_file
    tuple val(sample_id), path(first_read), path(second_read), emit: trimmed_fastqs
    path("${sample_id}_mapping_rate.csv"), emit: map_rate_ch

    script:
    """
    bowtie2 -p $threads -x $btidx -1 $first_read -2 $second_read | samtools sort -@ $threads -o ${sample_id}".sorted.bam"
    mapping_rate=\$(grep "overall alignment rate" .command.err | awk '{ print \$1 }')
    echo ${sample_id},\${mapping_rate} > ${sample_id}_mapping_rate.csv
    """
}

process get_overall_mapping_rate {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: 'mapping_rates.csv'
    input:
    file(mapping_rate)

    output:
    path("mapping_rates.csv")

    script:
    """
    echo "Sample,Mapping_Rate" > mapping_rates.csv
    cat *_mapping_rate.csv >> mapping_rates.csv
    """
}

process instrain {
    container '/software/pathogen/images/instrain-1.6.4-c2.simg'
    if (params.instrain_full_output) { publishDir path: "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*_instrain_output" }
    if (params.instrain_quick_profile) { publishDir path: "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*_instrain_quick_profile_output" }
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: '*.tsv'
    input:
    tuple val(sample_id), file(sorted_bam)
    val genome_file
    val stb_file
    val threads

    output:
    path ("${sample_id}_instrain_output"), optional: true
    path("${sample_id}_instrain_quick_profile_output"), optional: true
    path("${genome_info_file}"), optional: true
    path(sorted_bam), emit: sorted_bam
    path("${workdir}"), emit: workdir
    val sample_id, emit: sample_id

    script:
    genome_info_file="${sample_id}_genome_info.tsv"
    workdir="workdir.txt"
    """
    pwd > workdir.txt
    if $params.instrain_quick_profile
    then
        inStrain quick_profile $sorted_bam $genome_file -o ${sample_id}_instrain_quick_profile_output -p $threads -s $stb_file
    else
        inStrain profile $sorted_bam $genome_file -o ${sample_id}_instrain_output -p $threads -s $stb_file --database_mode --skip_plot_generation
    fi
    if ! $params.instrain_full_output && ! $params.instrain_quick_profile
    then
        mv ${sample_id}_instrain_output/output/${sample_id}"_instrain_output_genome_info.tsv" ./${sample_id}"_genome_info.tsv"
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
    get_overall_mapping_rate(bowtie2samtools.out.map_rate_ch.collect())
    if (!params.bowtie2_samtools_only) {
        instrain(bowtie2samtools.out.bam_file, params.genome_file, params.stb_file, params.instrain_threads)
    }
    if (!params.keep_bowtie2samtools && params.bowtie2_samtools_only) {
        cleanup_sorted_bam_files(bowtie2samtools.out.bam_file)
    }
    else if(!params.keep_bowtie2samtools) {
        cleanup_sorted_bam_files(instrain.out.sorted_bam)
    }
    if (!params.keep_instrain && !params.bowtie2_samtools_only) {
        cleanup_instrain_output(instrain.out.workdir, instrain.out.sample_id)
    }
}
