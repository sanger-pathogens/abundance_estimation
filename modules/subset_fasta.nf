process SUBSET_SOURMASH_GENES {
    input:
    tuple val(sample_id), path(sourmash_genes)

    output:
    tuple val(sample_id), path("${sample_id}_sourmash_genome*"), emit: subset_genomes
    script:
    """
    # select unique species from sourmash gather output
    cat ${sourmash_genes} | awk '{ print \$2,\$3 }' | sort -u > species.txt
    # split up for parallelisation in subsetting gtdb fasta file
    split -l 1 species.txt ${sample_id}_sourmash_genome
    """
}

process SUBSET_GTDB {
    maxForks 20 // limit so the gtdb fasta file isn't hit 100s of times at once
    input:
    tuple val(sample_id), path(sourmash_gene)

    output:
    tuple val(sample_id), path("*.fa"), emit: subset_fasta
    script:
    """
    /data/pam/team230/ol6/scratch/team162_pipelines/sourmash_test/abundance_estimation/splitter.sh $sourmash_gene $sample_id
    """
}

process CONCATENATE_FASTAS {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: '*_gtdb_subset.fa'
    input:
    tuple val(sample_id), path(fastas)

    output:
    tuple val(sample_id), path(subset_fasta), emit: subset_reference

    script:
    subset_fasta="${sample_id}_gtdb_subset.fa"
    """
    cat *.fa > ${sample_id}_gtdb_subset.fa
    """
}