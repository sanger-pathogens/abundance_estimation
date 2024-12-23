process MERGE_FASTQS {
    tag "${sample_id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_long'

    container 'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0'
    
    input:
    tuple val(sample_id), path(first_read), path(second_read)

    output:
    tuple val(sample_id), path(merged_fastq), emit: merged_fastq

    script:
    merged_fastq="${sample_id}_merged.fastq.gz"
    """
    echo "START MERGING"
    ls -l ${first_read}
    ls -l ${second_read}
    cat ${first_read} ${second_read} > ${merged_fastq}
    echo "FINISHED MERGING"
    ls -l ${merged_fastq}
    """
}
