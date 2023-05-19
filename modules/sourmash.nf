process SOURMASH_SKETCH {
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
    input:
    path(merged_fastq)

    output:
    path(sourmash_sketch), emit: sketch

    script:
    sourmash_sketch="SOURMASH_SKETCH"
    """
    sourmash sketch dna -p scaled=10000,k=31 ${merged_fastq} -o SOURMASH_SKETCH
    """
}

process SOURMASH_GATHER {
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
    input:
    path(sourmash_sketch)

    output:
    path(sourmash_genomes), emit: sourmash_genomes

    script:
    sourmash_genomes="sourmash_genomes.txt"
    """
    sourmash gather --dna ${sourmash_sketch} /data/pam/team162/shared/sourmash_db/gtdb-rs207.genomic-reps.dna.k31.zip -o sourmash.out
    # get species names out of sourmash output
    tail -n +2 sourmash.out | awk -F "," '{ print \$10 }' | sed 's|[][]||g' | sed 's|"||g' | awk '{ print \$1 }' > sourmash_genomes.txt
    """
}