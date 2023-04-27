process GET_GENOME_SUBSET {
    input:
    tuple val(sample_id), path(sourmash_genes)

    output:
    tuple val(sample_id), path("${sample_id}_ids.txt"), emit: subset_genomes
    script:
    """
    # select unique ids from sourmash gather output
    grep -f ${sourmash_genes} /lustre/scratch125/pam/teams/team230/ol6/team162_pipelines/sourmash_test/abundance_estimation/gtdb_r207_genomes.txt | sort -u > ${sample_id}_ids.txt
    """
}

process SUBSET_GTDB {
    publishDir "${params.results_dir}/subset_fastas/${sample_id}", mode: 'copy', overwrite: true, pattern: '*_gtdb_subset.fa'
    container '/software/pathogen/images/seqkit-2.3.1--h9ee0642_0.simg'
    maxForks 20
    input:
    tuple val(sample_id), path(sourmash_gene)

    output:
    tuple val(sample_id), path(subset_fasta), emit: subset_fasta
    script:
    subset_fasta="${sample_id}_gtdb_subset.fa"
    """
    seqkit grep -n -f ${sample_id}_ids.txt /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.fasta > ${sample_id}_gtdb_subset.fa
    """
}