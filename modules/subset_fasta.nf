process SUBSET_GENOMES {
    label 'cpu_1'
    label 'mem_50M'
    label 'time_queue_from_normal'

    input:
    tuple val(sample_id), path(sourmash_genomes)

    output:
    tuple val(sample_id), path("subset_ref_database.fasta"), emit: subset_genome

    script:
    genome_dir=params.genome_dir
    genomes_file_ext=params.genomes_file_ext
    aws_cli=params.aws_cli
    """
    subset_fasta.py ${genome_dir} ${sourmash_genomes} 
    """
}