process SOURMASH_SKETCH {
    tag "${sample_id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'

    container 'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0'

    input:
    tuple val(sample_id), path(merged_fastq)

    output:
    tuple val(sample_id), path(sourmash_sketch), emit: sketch

    script:
    sourmash_sketch="${sample_id}_sourmash_sketch"
    """
    sourmash sketch dna -p scaled=10000,k=31 ${merged_fastq} -o ${sample_id}_sourmash_sketch
    """
}

process SOURMASH_GATHER {
    tag "${sample_id}"
    label 'cpu_1'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0'

    input:
    tuple val(sample_id), path(sourmash_sketch)

    output:
    tuple val(sample_id), path(sourmash_genomes), emit: sourmash_genomes

    script:
    sourmash_genomes="${sample_id}_sourmash_genomes.txt"
    """
    sourmash gather --dna ${sourmash_sketch} ${params.sourmash_db} -o sourmash.out
    # get species names out of sourmash output
    tail -n +2 sourmash.out | awk -F "," '{ print \$10 }' | sed 's|[][]||g' | sed 's|"||g' | awk '{ print \$1 }' > ${sample_id}_sourmash_genomes.txt
    """
}