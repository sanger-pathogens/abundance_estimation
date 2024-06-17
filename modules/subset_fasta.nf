process SUBSET_GENOMES {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'

    input:
    tuple val(sample_id), path(sourmash_genomes)

    output:
    tuple val(sample_id), path("subset_ref_database.fasta"), emit: subset_genome
    script:
    """
    ${projectDir}/bin/subset_fasta.sh ${sourmash_genomes} ${params.genome_dir} ${params.genomes_file_ext}
    """
}