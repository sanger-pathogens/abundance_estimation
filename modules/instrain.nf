process INSTRAIN {
    tag "${sample_id}"
    container '/software/pathogen/images/instrain-1.6.4-c2.simg'

    input:
    tuple val(sample_id), path(sorted_bam), path(stb_file), path(genome_file)

    output:
    tuple val(sample_id), path("*_genome_info.tsv"), emit: genome_info_ch

    script:
    """
    assembly=\$(ls *.bam | sed 's/^.*GCF/GCF/' | sed 's/^.*GCA/GCA/' | sed 's/.sorted.bam//g')
    # run inStrain
    inStrain profile $sorted_bam $genome_file -o ${sample_id}_\${assembly}_instrain_output -p ${params.instrain_threads} -s $stb_file --database_mode --skip_plot_generation
    tail -n +2 ${sample_id}_\${assembly}_instrain_output/output/${sample_id}_\${assembly}_instrain_output_genome_info.tsv > ${sample_id}_\${assembly}_genome_info.tsv
    """
}

process GENERATE_STB {
    tag "${sample_id}"
    container '/software/pathogen/images/drep-3.2.2-c2.simg'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path(sorted_bam), path("*.stb"), path("*.fa"), emit: instrain_ch

    script:
    """
    assembly=\$(ls *.bam | sed 's/^.*GCF/GCF/' | sed 's/^.*GCA/GCA/' | sed 's/.sorted.bam//g')
    gunzip -c /data/pam/team162/shared/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207_genome_dir/\${assembly}_genomic.fna.gz > \${assembly}.fa
    parse_stb.py --reverse -f \${assembly}.fa -o \${assembly}_representative.stb
    """
}

process COLLATE_INSTRAIN_RESULTS {
    tag "${sample_id}"
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: '${sample_id}_genome_info.tsv'
    input:
    tuple val(sample_id), path(genome_info_files)
    path(instrain_header_file)

    output:
    tuple val(sample_id), path(instrain_results_file), emit: results_ch
    script:
    instrain_results_file="${sample_id}_genome_info.tsv"
    """
    cat $instrain_header_file > instrain_results.tsv
    cat *_genome_info.tsv >> instrain_results.tsv
    mv instrain_results.tsv ${sample_id}_genome_info.tsv
    """
}