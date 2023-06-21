process SUBSET_GTDB {
    input:
    tuple val(sample_id), path(sourmash_genomes)

    output:
    tuple val(sample_id), path("subset_gtdb_ref.fasta"), emit: subset_genome
    script:
    """
    while read genome
    do
      zcat ${params.genome_dir}/\${genome}_genomic.fna.gz >> subset_gtdb_ref.fasta
    done < ${sourmash_genomes}
    """
}