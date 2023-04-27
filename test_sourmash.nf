#!/usr/bin/env nextflow

include { MERGE_FASTQS } from './modules/merge_fastq.nf'
include { SOURMASH_SKETCH; SOURMASH_GATHER } from './modules/sourmash.nf'
include { SUBSET_SOURMASH_GENES; SUBSET_GTDB; CONCATENATE_FASTAS } from './modules/subset_fasta.nf'
include { BOWTIE_INDEX } from './modules/bowtie.nf'

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.first_read), file(row.second_read)) }

    MERGE_FASTQS(fastq_path_ch)

    SOURMASH_SKETCH(MERGE_FASTQS.out.merged_fastq)

    SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch)

    SUBSET_SOURMASH_GENES(SOURMASH_GATHER.out.sourmash_genomes)

    SUBSET_GTDB(SUBSET_SOURMASH_GENES.out.subset_genomes.transpose())

    CONCATENATE_FASTAS(SUBSET_GTDB.out.subset_fasta.groupTuple())

    BOWTIE_INDEX(CONCATENATE_FASTAS.out.subset_reference)
}
