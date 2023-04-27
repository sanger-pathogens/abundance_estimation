#!/usr/bin/env nextflow

include { MERGE_FASTQS } from './modules/merge_fastq.nf'
include { SOURMASH_SKETCH; SOURMASH_GATHER } from './modules/sourmash.nf'
include { GET_GENOME_SUBSET; SUBSET_GTDB } from './modules/subset_fasta.nf'
include { BOWTIE_INDEX } from './modules/bowtie.nf'

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.first_read), file(row.second_read)) }

    MERGE_FASTQS(fastq_path_ch)

    SOURMASH_SKETCH(MERGE_FASTQS.out.merged_fastq)

    SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch)

    GET_GENOME_SUBSET(SOURMASH_GATHER.out.sourmash_genomes)

    SUBSET_GTDB(GET_GENOME_SUBSET.out.subset_genomes)

    BOWTIE_INDEX(SUBSET_GTDB.out.subset_fasta)

    //BOWTIE_INDEX.out.bowtie_index
    //                           .join(SUBSET_GTDB.out.subset_fasta)
    //                           .join(fastq_path_ch)
}
