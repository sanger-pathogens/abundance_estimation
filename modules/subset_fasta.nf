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
    while read genome
    do
      zcat -f ${params.genome_dir}/\${genome}${params.genomes_file_ext} >> subset_ref_database.fasta
    done < ${sourmash_genomes}
    """
}