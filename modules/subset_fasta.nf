process SUBSET_GENOMES {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'

    input:
    tuple val(sample_id), path(sourmash_genomes)

    output:
    tuple val(sample_id), path("subset_ref_database.fasta"), emit: subset_genome

    script:
    genome_dir=params.genome_dir
    genomes_file_ext=params.genomes_file_ext
    aws_cli=params.aws_cli
    template "subset_fasta.sh"
}