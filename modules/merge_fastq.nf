process MERGE_FASTQS {
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
    input:
    path(fastqs)

    output:
    path(merged_fastq), emit: merged_fastq

    script:
    merged_fastq="merged.fastq.gz"
    """
    cat *.fastq.gz > merged.fastq.gz
    """
}
