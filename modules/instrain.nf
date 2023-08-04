process INSTRAIN {
    tag "${sample_id}"
    label 'time_queue_from_normal'

    container '/software/pathogen/images/instrain-1.6.4-c2.simg'

    if (params.instrain_full_output) { publishDir path: "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*_instrain_output" }
    if (params.instrain_quick_profile) { publishDir path: "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*_instrain_quick_profile_output" }
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: '*.tsv'
    
    input:
    tuple val(sample_id), path(sorted_bam), path(stb_file), path(genome_file)

    output:
    path ("${sample_id}_instrain_output"), optional: true
    path("${sample_id}_instrain_quick_profile_output"), optional: true
    path("${genome_info_file}"), optional: true
    path(sorted_bam), emit: sorted_bam
    path("${workdir}"), emit: workdir
    val sample_id, emit: sample_id

    script:
    genome_info_file="${sample_id}_genome_info.tsv"
    workdir="workdir.txt"
    """
    pwd > workdir.txt
    if $params.instrain_quick_profile
    then
        inStrain quick_profile $sorted_bam $genome_file -o ${sample_id}_instrain_quick_profile_output -p ${task.cpus} -s $stb_file
    else
        inStrain profile $sorted_bam $genome_file -o ${sample_id}_instrain_output -p ${task.cpus} -s $stb_file --database_mode --skip_plot_generation
    fi
    if ! $params.instrain_full_output && ! $params.instrain_quick_profile
    then
        mv ${sample_id}_instrain_output/output/${sample_id}"_instrain_output_genome_info.tsv" ./${sample_id}"_genome_info.tsv"
    fi
    """
}

process GENERATE_STB {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'

    input:
    tuple val(sample_id), path(sourmash_genomes)

    output:
    tuple val(sample_id), path("*.stb"), emit: stb_ch

    script:
    """
    sed 's|\$|_genomic.fna.gz|g' $sourmash_genomes > genomes.txt
    grep -w -f genomes.txt ${params.stb_file} > ${sample_id}_gtdb_subset.stb
    """
}
