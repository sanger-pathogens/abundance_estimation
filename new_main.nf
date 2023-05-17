#!/usr/bin/env nextflow

include { METAWRAP_QC } from './modules/metawrap_qc.nf'
include { CLEANUP_SORTED_BAM_FILES; CLEANUP_TRIMMED_FASTQ_FILES; CLEANUP_INSTRAIN_OUTPUT } from './modules/cleanup.nf'
include { VALIDATE_PARAMETERS; PRINT_HELP } from './modules/helper_functions.nf'
include { MERGE_FASTQS } from './modules/merge_fastq.nf'
include { SOURMASH_SKETCH; SOURMASH_GATHER; SORT_SOURMASH_GENOMES } from './modules/sourmash.nf'
include { SUBSET_GTDB } from './modules/subset_fasta.nf'
include { BOWTIE_INDEX; BOWTIE2SAMTOOLS; GET_OVERALL_MAPPING_RATE } from './modules/bowtie.nf'
include { GENERATE_STB; INSTRAIN } from './modules/instrain.nf'

workflow {
    if (params.help) {
        PRINT_HELP()
        exit 0
    }

    VALIDATE_PARAMETERS()

    manifest_ch = Channel.fromPath(params.manifest)

    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.first_read), file(row.second_read)) }

    MERGE_FASTQS(fastq_path_ch)

    SOURMASH_SKETCH(MERGE_FASTQS.out.merged_fastq)

    SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch)

    SORT_SOURMASH_GENOMES(SOURMASH_GATHER.out.sourmash_genomes.collect())

    SUBSET_GTDB(SORT_SOURMASH_GENOMES.out.sorted_genomes)

    BOWTIE_INDEX(SUBSET_GTDB.out.subset_genome)

    GENERATE_STB(SORT_SOURMASH_GENOMES.out.sorted_genomes)

    if (params.skip_qc) {
        BOWTIE2SAMTOOLS(fastq_path_ch, BOWTIE_INDEX.out.bowtie_index, params.bowtie2_samtools_threads)
    }
    else {
        METAWRAP_QC(fastq_path_ch)
        BOWTIE2SAMTOOLS(METAWRAP_QC.out.trimmed_fastqs, BOWTIE_INDEX.out.bowtie_index, params.bowtie2_samtools_threads)
    }
    if (!params.keep_metawrap_qc) {
        if (!params.skip_qc) {
            CLEANUP_TRIMMED_FASTQ_FILES(BOWTIE2SAMTOOLS.out.trimmed_fastqs)
        }
    }
    GET_OVERALL_MAPPING_RATE(BOWTIE2SAMTOOLS.out.map_rate_ch.collect())
    if (!params.bowtie2_samtools_only) {
        INSTRAIN(BOWTIE2SAMTOOLS.out.bam_file, SUBSET_GTDB.out.subset_genome, GENERATE_STB.out.stb_ch, params.instrain_threads)
    }
    if (!params.keep_bowtie2samtools && params.bowtie2_samtools_only) {
        CLEANUP_SORTED_BAM_FILES(BOWTIE2SAMTOOLS.out.bam_file)
    }
    else if(!params.keep_bowtie2samtools) {
        CLEANUP_SORTED_BAM_FILES(INSTRAIN.out.sorted_bam)
    }
    if (!params.keep_instrain && !params.bowtie2_samtools_only) {
        CLEANUP_INSTRAIN_OUTPUT(INSTRAIN.out.workdir, INSTRAIN.out.sample_id)
    }
}