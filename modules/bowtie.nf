process BOWTIE_INDEX {
    container '/software/pathogen/images/bowtie2-2.3.5--py37he860b03_0.simg'
    publishDir "${params.results_dir}/bowtie_indexes/${sample_id}", mode: 'copy', overwrite: true, pattern: '*bt2*'
    input:
    tuple val(sample_id), path(subset_fasta)

    output:
    tuple val(sample_id), path(subset_fasta), path("*.bt2"), emit: bowtie_index

    script:
    """
    bowtie2-build $subset_fasta ${sample_id}.fasta --threads 16 --large-index
    """
}