#!/usr/bin/env nextflow

include { MERGE_FASTQS } from './modules/merge_fastq.nf'
include { SOURMASH_SKETCH; SOURMASH_GATHER } from './modules/sourmash.nf'
include { SUBSET_GTDB } from './modules/subset_fasta.nf'
include { BOWTIE_INDEX } from './modules/bowtie.nf'

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.first_read), file(row.second_read)) }

    fastq_ch = fastq_path_ch.map{ sample_id, first_read, second_read -> tuple(file(first_read), file(second_read)) }

    //fastq_ch.collect().view()
    MERGE_FASTQS(fastq_ch.collect())

    SOURMASH_SKETCH(MERGE_FASTQS.out.merged_fastq)

    SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch)

    SUBSET_GTDB(SOURMASH_GATHER.out.sourmash_genomes)

    BOWTIE_INDEX(SUBSET_GTDB.out.subset_genome)

    //BOWTIE_INDEX.out.bowtie_index
    //                           .join(SUBSET_GTDB.out.subset_fasta)
    //                           .join(fastq_path_ch)
}
