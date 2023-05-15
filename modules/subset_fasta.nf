process SUBSET_ASSEMBLY_FILE {
    tag "${sample_id}"
    input:
    tuple val(sample_id), path(sourmash_genomes_csv)

    output:
    path("${sample_id}_genomes_to_be_indexed.csv"), emit: assemblies_to_be_indexed
    path(sourmash_genomes_csv), emit: sourmash_genomes
    script:
    """
    # check what species are already indexed in cache
    ls ${params.index_cache}/*/*.bt2 | awk -F "." '{ print \$1 }' | sort -u | rev | cut -d"/" -f1 | rev > existing_genomes.txt
    grep -vf existing_genomes.txt ${sourmash_genomes_csv} > ${sample_id}_genomes_to_be_indexed.csv
    """
}

process GUNZIP_ASSEMBLY {
    tag "${sourmash_genome}"
    input:
    val(sourmash_genome)

    output:
    tuple val(sourmash_genome), path("*.fna"), emit: assembly

    script:
    """
    gunzip -c /data/pam/team162/shared/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207_genome_dir/${sourmash_genome}_genomic.fna.gz > ${sourmash_genome}.fna
    """
}

//process MERGE_FASTAS {
//    tag "${sample_id}"
//    input:
//    tuple val(sample_id), path(bams)
//
//    output:
//    tuple val(sample_id), path("${sample_id}_merged.fa")
//
//    script:
//    """
//    assemblies=\$(ls *.bam | sed 's/^.*GCF/GCF/' | sed 's/.sorted.bam//g')
//    for assembly in \${assemblies}
//    do
//      zcat /data/pam/team162/shared/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207_genome_dir/\${assembly}_genomic.fna.gz >> ${sample_id}_merged.fa
//    done
//    """
//}


//process GET_GENOME_SUBSET {
//    tag "${sourmash_genome}"
//    input:
//    val(sourmash_genome)
//
//    output:
//    tuple path("${sourmash_genome}_gtdb_ids.txt"), val(sourmash_genome), emit: gtdb_genome
//    script:
//    """
//    grep -h ${sourmash_genome} /lustre/scratch125/pam/teams/team230/ol6/team162_pipelines/sourmash_test/abundance_estimation/gtdb_r207_genomes.txt | sort -u > ${sourmash_genome}_gtdb_ids.txt
//    """
//}
//
//process SUBSET_GTDB {
//    tag "${sourmash_species}"
//    container '/software/pathogen/images/seqkit-2.3.1--h9ee0642_0.simg'
//    maxForks 20
//    input:
//    tuple path(genome_subset), val(sourmash_species)
//
//    output:
//    tuple val(sourmash_species), path(subset_fasta), emit: subset_fasta
//    script:
//    subset_fasta="${sourmash_species}_gtdb_subset.fa"
//    """
//    seqkit grep -n -f ${genome_subset} /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.fasta > ${sourmash_species}_gtdb_subset.fa
//    """
//}