#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { metawrap_qc } from './modules/metawrap_qc.nf'

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
    def output_list = ["0", "1"]
    if (params.full_output) {
        if (!output_list.any { it.contains(params.full_output.toString()) }) {
            log.error("Please specify the output as 0 or 1 using the --type option.")
            errors += 1
        }
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

    script:
    """
    bowtie2 -p $threads -x $btidx -1 $first_read -2 $second_read | samtools sort -@ $threads -o ${sample_id}".sorted.bam"
    """
}

process instrain {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true
    input:
    tuple val(sample_id), file(sorted_bam)
    file(genome_file)
    file(stb_file)
    val threads

    output:
    path("${genome_info_file}")

    script:
    genome_info_file="${sample_id}_genome_info.tsv"
    """
    inStrain profile $sorted_bam $genome_file -o $sample_id -p $threads -s $stb_file --database_mode --skip_plot_generation
    mv ${sample_id}/output/${sample_id}"_genome_info.tsv" .
    """
}

process instrain_full_output {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true
    input:
    tuple val(sample_id), file(sorted_bam)
    file(genome_file)
    file(stb_file)
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
    genome_file = Channel.fromPath(params.genome_file)
    stb_file = Channel.fromPath(params.stb_file)
    metawrap_qc(fastq_path_ch)
    bowtie2samtools(metawrap_qc.out.trimmed_fastqs, params.btidx, params.bowtie2_samtools_threads)
    if (params.full_output.toString() == "1") {
        instrain_full_output(bowtie2samtools.out.bam_file, genome_file, stb_file, params.instrain_threads)
    }
    else {
        instrain(bowtie2samtools.out.bam_file, genome_file, stb_file, params.instrain_threads)
    }
}
