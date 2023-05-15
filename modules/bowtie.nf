process BOWTIE_INDEX {
    tag "${sourmash_genome}"
    container '/software/pathogen/images/bowtie2-2.3.5--py37he860b03_0.simg'
    publishDir "${params.index_cache}/${sourmash_genome}", mode: 'copy', overwrite: true, pattern: '*bt2*'
    input:
    tuple val(sourmash_genome), path(assembly)

    output:
    path("*.bt2*"), emit: bowtie_index
    val(sourmash_genome), emit: genome

    script:
    """
    bowtie2-build $assembly ${sourmash_genome}.bt2 --threads 4
    """
}

//process DUMMY {
//    tag "${sample_id}"
//    input:
//    tuple val(sample_id), path(read_1), path(read_2), val(assembly_id)
//    path(indexes)
//
//    output:
//
//    script:
//    """
//    echo ${sample_id} ${read_1} ${read_2} ${assembly_id}
//    echo ${indexes}
//    """
//}

process BOWTIE2SAMTOOLS {
    tag "${sample_id}"
    container '/software/pathogen/images/bowtie2-samtools-1.0.simg'
    input:
    tuple val(sample_id), path(read_1), path(read_2), val(assembly_id)
    path(indexes)

    output:
    tuple val(sample_id), path("${sample_id}_${assembly_id}.sorted.bam"), emit: bam_file

    script:
    """
    bowtie2 -p ${params.bowtie2_samtools_threads} -x ${params.index_cache}/${assembly_id}/${assembly_id}.bt2 -1 $read_1 -2 $read_2 | samtools sort -@ ${params.bowtie2_samtools_threads} -o ${sample_id}_${assembly_id}".sorted.bam"
    """
}

//process SAMTOOLS_MERGE {
//    container '/software/pathogen/images/bowtie2-samtools-1.0.simg'
//    input:
//    tuple val(sample_id), path(bams)
//
//    output:
//    tuple val(sample_id), path("${sample_id}_merged.sorted.bam"), emit: bam_file
//
//    script:
//    """
//    samtools merge ${sample_id}_merged.sorted.bam -c -p *.bam
//    """
//}