def PRINT_HELP() {
    log.info """
    Usage:
    nextflow run .

    Options:
      --manifest                   Manifest containing paths to fastq files (mandatory)
      --bowtie2_samtools_threads   threads - default: 4 (optional)
      --instrain_threads           threads - default: 4 (optional)
      --instrain_full_output       get full instrain output - default false (optional)
      --cleanup_intermediate_files cleanup intermediate files - default false (optional)
      --skip_qc                    skip metawrap qc step - default false (optional)
      --stb_file                   stb file - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.stb (optional)
      --genome_dir                 genome folder - default: /data/pam/team162/shared/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207_genome_dir (optional)
      --results_dir                results folder - default: ./nextflow_results
      --sourmash_db                sourmash database - default: /data/pam/team162/shared/sourmash_db/gtdb-rs207.genomic-reps.dna.k31.zip (optional)
      --instrain_queue             job queue for instrain - default normal (optional)
      --instrain_quick_profile     use quick-profile option for inStrain - default false (optional)
      --bowtie2_samtools_only      only run bowtie2_samtools process - default false (optional)
      --help                       print this help message (optional)
    """.stripIndent()
}

def validate_path_param(
    param_option,
    param,
    type="file",
    mandatory=true) {
        valid_types=["file", "directory"]
        if (!valid_types.any { it == type }) {
                log.error("Invalid type '${type}'. Possibilities are ${valid_types}.")
                return 1
        }
        param_name = (param_option - "--").replaceAll("_", " ")
        if (param) {
            def file_param = file(param)
            if (!file_param.exists()) {
                log.error("The given ${param_name} '${param}' does not exist.")
                return 1
            } else if (
                (type == "file" && !file_param.isFile())
                ||
                (type == "directory" && !file_param.isDirectory())
            ) {
                log.error("The given ${param_name} '${param}' is not a ${type}.")
                return 1
            }
        } else if (mandatory) {
            log.error("No ${param_name} specified. Please specify one using the ${param_option} option.")
            return 1
        }
        return 0
    }

def validate_number_param(param_option, param) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param != null) /* Explicit comparison with null, because 0 is an acceptable value */ {
        if (!(param instanceof Number)) {
            log.error("The ${param_name} specified with the ${param_option} option must be a valid number")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

def validate_results_dir(results_dir) {
    results_dir = file(results_dir)
    if (results_dir.exists() && !results_dir.isDirectory()) {
        log.error("The given results_dir '${results_dir}' is not a directory.")
        return 1
    }
    return 0
}

def VALIDATE_PARAMETERS() {
    def errors = 0

    errors += validate_path_param("--manifest", params.manifest)
    errors += validate_path_param("--stb_file", params.stb_file)
    errors += validate_path_param("--sourmash_db", params.sourmash_db)
    errors += validate_path_param("--genome_dir", params.genome_dir, type="directory")
    errors += validate_number_param("--bowtie2_samtools_threads", params.bowtie2_samtools_threads)
    errors += validate_number_param("--queue_size", params.queue_size)
    errors += validate_number_param("--instrain_threads", params.instrain_threads)
    errors += validate_results_dir(params.results_dir)

    if (params.instrain_full_output && params.instrain_quick_profile) {
        log.error("the --instrain_full_output and --instrain_quick_profile options are incompatible, please choose one of these options.")
        errors += 1
    }

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}