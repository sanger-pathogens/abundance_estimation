#!/usr/bin/env nextflow

/*
========================================================================================
    HELP
========================================================================================
*/

def printHelp() {
    log.info """
    Usage:
    nextflow run main.nf

    Options:
      --manifest                      Manifest containing paths to fastq files. with headers ID,R1,R2. (mandatory)
      --outdir                        Name of results folder. [default: ./results] (optional)
      --bowtie2_samtools_threads      Threads for bowtie2 and samtools. [default: 4] (optional)
      --instrain_threads              Threads for instrain. [default: 4] (optional)
      --instrain_full_output          Get full instrain output. [default: false] (optional)
      --cleanup_intermediate_files    Cleanup intermediate files. [default: false] (optional)
      --skip_qc                       Skip metawrap qc step. [default: false] (optional)
      --stb_file                      Supply stb file. [default: /data/pam/software/inStrain/stb/gtdb_genomes_reps_r220.stb] (optional)
      --genome_dir                    Supply genome folder. [default: /data/pam/software/inStrain/genomes/gtdb_genomes_reps_r220] (optional)
      --sourmash_db                   Supply sourmash database. [default: /data/pam/software/sourmash/signatures_zipped/gtdb_genomes_reps_r220.zip] (optional)
      --instrain_quick_profile        Use quick-profile option for inStrain. [default: false] (optional)
      --bowtie2_samtools_only         Only run bowtie2_samtools process. [default: false] (optional)
      --help                          Print this help message. (optional)
    """.stripIndent()
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

    stb_ch = Channel.fromPath(params.stb_file)

    SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch,stb_ch)

    SUBSET_GENOMES(SOURMASH_GATHER.out.sourmash_genomes)

    BOWTIE_INDEX(SUBSET_GENOMES.out.subset_genome)

    SUBSET_STB(SOURMASH_GATHER.out.sourmash_genomes,stb_ch)

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
