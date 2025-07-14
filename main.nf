#!/usr/bin/env nextflow

/*
========================================================================================
    HELP
========================================================================================
*/

def printHelp() {
    NextflowTool.help_message("${workflow.ProjectDir}/schema.json", [],
    params.monochrome_logs, log)
}


if (params.help) {
    printHelp()
    exit 0
}

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { validate_parameters } from './modules/helper_functions.nf'
include { CLEANUP_SORTED_BAM_FILES; CLEANUP_TRIMMED_FASTQ_FILES; CLEANUP_INSTRAIN_OUTPUT } from './modules/cleanup.nf'
include { MERGE_FASTQS } from './modules/merge_fastq.nf'
include { SOURMASH_SKETCH; SOURMASH_GATHER } from './modules/sourmash.nf'
include { SUBSET_GENOMES } from './modules/subset_fasta.nf'
include { BOWTIE_INDEX; BOWTIE2SAMTOOLS; GET_OVERALL_MAPPING_RATE } from './modules/bowtie.nf'
include { SUBSET_STB; INSTRAIN } from './modules/instrain.nf'

//
// SUBWORKFLOWS
//
include { METAWRAP_QC } from './subworkflows/metawrap_qc.nf'

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

validate_parameters()

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.ID, file(row.R1), file(row.R2)) }

    if (params.skip_qc) {
        MERGE_FASTQS(fastq_path_ch)
    } else {
        METAWRAP_QC(fastq_path_ch)

        MERGE_FASTQS(METAWRAP_QC.out.filtered_reads)
    }

    SOURMASH_SKETCH(MERGE_FASTQS.out.merged_fastq)

    /*stb_ch = Channel.fromPath(params.stb_file)*/

    SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch,file(params.stb_file)/*stb_ch*/)

    SUBSET_GENOMES(SOURMASH_GATHER.out.sourmash_genomes)

    BOWTIE_INDEX(SUBSET_GENOMES.out.subset_genome)

    SUBSET_STB(SOURMASH_GATHER.out.sourmash_genomes,file(params.stb_file)/*stb_ch*/)

    if (params.skip_qc) {
        bowtie_mapping_ch = fastq_path_ch.join(BOWTIE_INDEX.out.bowtie_index)
        BOWTIE2SAMTOOLS(bowtie_mapping_ch, params.bowtie2_samtools_threads)
    }
    else {
        bowtie_mapping_ch = METAWRAP_QC.out.filtered_reads.join(BOWTIE_INDEX.out.bowtie_index)
        BOWTIE2SAMTOOLS(bowtie_mapping_ch, params.bowtie2_samtools_threads)
    }

    if (params.cleanup_intermediate_files) {
        if (!params.skip_qc) {
            CLEANUP_TRIMMED_FASTQ_FILES(BOWTIE2SAMTOOLS.out.trimmed_fastqs)
        }
    }

    GET_OVERALL_MAPPING_RATE(BOWTIE2SAMTOOLS.out.map_rate_ch.collect())
 
    if (!params.bowtie2_samtools_only) {
        instrain_profiling_ch = BOWTIE2SAMTOOLS.out.bam_file.join(SUBSET_STB.out.sub_stb_ch).join(SUBSET_GENOMES.out.subset_genome)
        INSTRAIN(instrain_profiling_ch)
    }

    if (params.cleanup_intermediate_files && params.bowtie2_samtools_only) {
        bam_files = BOWTIE2SAMTOOLS.out.bam_file.map { it[1] }
        CLEANUP_SORTED_BAM_FILES(bam_files)
    }
    else if(params.cleanup_intermediate_files) {
        CLEANUP_SORTED_BAM_FILES(INSTRAIN.out.sorted_bam)
    }

    if (params.cleanup_intermediate_files && !params.bowtie2_samtools_only) {
        CLEANUP_INSTRAIN_OUTPUT(INSTRAIN.out.workdir, INSTRAIN.out.sample_id)
    }
}
