process SOURMASH_SKETCH {
    tag "${sample_id}"
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
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
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
    input:
    tuple val(sample_id), path(sourmash_sketch)

    output:
    path(sourmash_genomes), emit: sourmash_genomes

    script:
    sourmash_genomes="${sample_id}_sourmash_genomes.txt"
    """
    sourmash gather --dna ${sourmash_sketch} /data/pam/team162/shared/sourmash_db/gtdb-rs207.genomic-reps.dna.k31.zip -o sourmash.out
    # get species names out of sourmash output
    tail -n +2 sourmash.out | awk -F "," '{ print \$10 }' | sed 's|[][]||g' | sed 's|"||g' | awk '{ print \$1 }' > ${sample_id}_sourmash_genomes.txt
    """
}

process SORT_SOURMASH_GENOMES {
    input:
    path(sourmash_genomes)

    output:
    path(sorted_sourmash_genomes), emit: sorted_genomes

    script:
    sorted_sourmash_genomes="sorted_sourmash_genomes.txt"
    """
    cat *.txt | sort -u > sorted_sourmash_genomes.txt
    """
}