#!/usr/bin/env nextflow

include { MERGE_FASTQS } from './modules/merge_fastq.nf'
include { SOURMASH_SKETCH; SOURMASH_GATHER } from './modules/sourmash.nf'
include { SUBSET_ASSEMBLY_FILE; GUNZIP_ASSEMBLY; COLLECT_ASSEMBLIES } from './modules/subset_fasta.nf'
include { BOWTIE_INDEX; BOWTIE2SAMTOOLS } from './modules/bowtie.nf'
include { INSTRAIN; GENERATE_STB; COLLATE_INSTRAIN_RESULTS } from './modules/instrain.nf'

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.first_read), file(row.second_read)) }

    MERGE_FASTQS(fastq_path_ch)

    SOURMASH_SKETCH(MERGE_FASTQS.out.merged_fastq)

    SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch)

    SUBSET_ASSEMBLY_FILE(SOURMASH_GATHER.out.sourmash_genomes)

    sourmash_assembly_channel = SUBSET_ASSEMBLY_FILE.out.sourmash_genomes.splitCsv(header: true, sep: ',')
        .map{ row -> tuple(row.sample_id, row.genome) }
        .groupTuple()

    indexing_ch = SUBSET_ASSEMBLY_FILE.out.assemblies_to_be_indexed.splitCsv(header: true, sep: ',')
        .map{ row -> tuple(row.genome) }
        .collect()
        .flatten()
        .unique()

    bowtie_index_ch = Channel.empty()

    COLLECT_ASSEMBLIES(indexing_ch)

    BOWTIE_INDEX(COLLECT_ASSEMBLIES.out.assembly)

    //bowtie_index_ch = bowtie_index_ch.mix(BOWTIE_INDEX.out.bowtie_index.last())

    // if no indexes are computed, add a file into the channel to avoid dependency issues and proceed straight to mapping
    bowtie_index_ch = BOWTIE_INDEX.out.bowtie_index.ifEmpty(file(params.index_cache))

    //mapping_ch = MERGE_FASTQS.out.fastq_ch.join(sourmash_assembly_channel).transpose()

    mapping_ch = sourmash_assembly_channel.transpose()

    BOWTIE2SAMTOOLS(mapping_ch, bowtie_index_ch.last(), file(params.manifest))

    GENERATE_STB(BOWTIE2SAMTOOLS.out.bam_file)

    INSTRAIN(GENERATE_STB.out.instrain_ch)

    instrain_genome_info_ch = INSTRAIN.out.genome_info_ch.groupTuple()

    COLLATE_INSTRAIN_RESULTS(instrain_genome_info_ch, file(params.instrain_header_file))
//
    //SAMTOOLS_MERGE(BOWTIE2SAMTOOLS.out.bam_file.groupTuple())

    //MERGE_FASTAS(BOWTIE2SAMTOOLS.out.bam_file.groupTuple())


    //GET_GENOME_SUBSET(indexing_ch)

    //SUBSET_GTDB(GET_GENOME_SUBSET.out.gtdb_species)

    //BOWTIE_INDEX(SUBSET_GTDB.out.subset_fasta)


    // BOWTIE2 INDEX
    //ref_without_extension = "${reference.parent}/${reference.baseName}"
    //bt2_index_files = file("${ref_without_extension}*.bt2")
    //if (bt2_index_files) {
    //    Channel.fromPath(bt2_index_files)
    //        .collect()
    //        .map { bt2_index_files -> tuple(ref_without_extension, bt2_index_files) }
    //        .dump(tag: 'bt2_index')
    //        .set { ch_bt2_index }
    //} else {
    //    BOWTIE2_INDEX(
    //        reference
    //    )
    //    BOWTIE2_INDEX.out.bt2_index.dump(tag: 'bt2_index').set { ch_bt2_index }
    //}


    //ENTREZ_TEST(SOURMASH_GATHER.out.sourmash_genomes)

    //GET_GENOME_SUBSET(SOURMASH_GATHER.out.sourmash_genomes)

    //GET_GENOME_SUBSET.out.subset_genomes.collect().view()

    //SUBSET_GTDB(GET_GENOME_SUBSET.out.subset_genomes)

    //BOWTIE_INDEX(SUBSET_GTDB.out.subset_fasta
}
