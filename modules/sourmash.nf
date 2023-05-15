process SOURMASH_SKETCH {
    tag "${sample_id}"
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
    tag "${sample_id}"
    container '/software/pathogen/images/sourmash-4.5.0--hdfd78af_0.simg'
    input:
    tuple val(sample_id), path(sourmash_sketch)

    output:
    tuple val(sample_id), path(sourmash_genomes), emit: sourmash_genomes

    script:
    sourmash_genomes="${sample_id}_sourmash_genomes.csv"
    """
    sourmash gather --dna ${sourmash_sketch} /data/pam/team162/shared/sourmash_db/gtdb-rs207.genomic-reps.dna.k31.zip -o sourmash.out
    # get genomes out of sourmash output
    tail -n +2 sourmash.out | head -3 | awk -F "," '{ print \$10 }' | awk '{ print \$1 }' | sed 's|"||g' > ${sample_id}_sourmash_genomes.txt
    echo "sample_id,genome" > ${sample_id}_sourmash_genomes.csv
    sed "s|^|${sample_id},|g" ${sample_id}_sourmash_genomes.txt >> ${sample_id}_sourmash_genomes.csv
    """
}

//process ENTREZ_TEST {
//    maxForks 1
//    container '/software/pathogen/images/entrez-direct-16.2--he881be0_1.simg'
//    input:
//    tuple val(sample_id), path(sourmash_genomes)
//
//    output:
//    tuple val(sample_id), path(subset_genomes), emit: subset_genomes
//
//    script:
//    subset_genomes="${sample_id}_subset_genomes.txt"
//    """
//    split -l 1 $sourmash_genomes sourmash
//    rm $sourmash_genomes
//    for f in sourmash*
//    do
//        pattern=\$(cat \$f)
//        esearch -db assembly -query \${pattern} | elink -db assembly -target nuccore | efetch | grep -A10 "Seq-entry" | grep accession | awk '{ print \$NF }' | sed 's|"||g' | sed 's|,||g' > tmp_query_file.txt
//        grep -f tmp_query_file.txt /lustre/scratch125/pam/teams/team230/ol6/team162_pipelines/sourmash_test/abundance_estimation/gtdb_r207_genomes.txt >> ${sample_id}_subset_genomes.txt
//        sleep 60
//    done
//    rm tmp_query_file.txt
//    """
//}