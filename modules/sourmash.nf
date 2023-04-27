process SOURMASH_SKETCH {
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
    input:
    tuple val(sample_id), path(merged_fastq)

    output:
    tuple val(sample_id), path(SOURMASH_SKETCH), emit: sketch

    script:
    SOURMASH_SKETCH="${sample_id}_SOURMASH_SKETCH"
    """
    sourmash sketch dna -p scaled=10000,k=31 ${merged_fastq} -o ${sample_id}_SOURMASH_SKETCH
    """
}

process SOURMASH_GATHER {
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
    input:
    tuple val(sample_id), path(sourmash_sketch)

    output:
    tuple val(sample_id), path(sourmash_genomes), emit: sourmash_genomes

    script:
    sourmash_genomes="${sample_id}_sourmash_genomes.txt"
    """
    sourmash gather --dna ${sourmash_sketch} /data/pam/team162/shared/sourmash_db/gtdb-rs207.genomic-reps.dna.k31.zip -o sourmash.out
    # this is arbritary for testing
    head sourmash.out | tail -n +2 sourmash.out | awk -F "," '{ print \$10 }' > ${sample_id}_sourmash_genomes.txt
    """
}