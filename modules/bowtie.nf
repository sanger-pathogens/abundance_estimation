process BOWTIE_INDEX {
    label 'cpu_4'
    label 'mem_32'
    label 'time_queue_from_normal'

    container '/software/pathogen/images/bowtie2-samtools-1.1-c1.simg'

    input:
    tuple val(sample_id), path(subset_fasta)

    output:
    tuple val(sample_id), path("*.bt2*"), emit: bowtie_index

    script:
    """
    bowtie2-build $subset_fasta ${sample_id}_gtdb_subset.bt2 --threads 4 --large-index
    """
}

process BOWTIE2SAMTOOLS {
    tag "${sample_id}"
    label 'time_queue_from_normal'

    container '/software/pathogen/images/bowtie2-samtools-1.1-c1.simg'

    if (params.bowtie2_samtools_only) { publishDir path: "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*.sorted.bam" }
    input:
    tuple val(sample_id), path(first_read), path(second_read), path(btidx)
    val threads

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam_file
    tuple val(sample_id), path(first_read), path(second_read), emit: trimmed_fastqs
    path("${sample_id}_mapping_rate.csv"), emit: map_rate_ch

    script:
    """
    bowtie2 -p $threads -x ${sample_id}_gtdb_subset.bt2 -1 $first_read -2 $second_read | samtools sort -@ $threads -o ${sample_id}".sorted.bam"
    mapping_rate=\$(grep "overall alignment rate" .command.err | awk '{ print \$1 }')
    echo ${sample_id},\${mapping_rate} > ${sample_id}_mapping_rate.csv
    """
}

process GET_OVERALL_MAPPING_RATE {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: 'mapping_rates.csv'
    
    input:
    path(mapping_rate)

    output:
    path("mapping_rates.csv")

    script:
    """
    echo "Sample,Mapping_Rate" > mapping_rates.csv
    cat *_mapping_rate.csv >> mapping_rates.csv
    """
}