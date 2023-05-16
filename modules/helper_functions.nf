def PRINT_HELP() {
    log.info """
    Usage:
    nextflow run .

    Options:
      --manifest                   Manifest containing paths to fastq files (mandatory)
      --bowtie2_samtools_threads   threads - default: 8 (optional)
      --instrain_threads           threads - default: 8 (optional)
      --instrain_full_output       get full instrain output - default false (optional)
      --keep_metawrap_qc           don't cleanup metawrap_qc output - default false (optional)
      --keep_bowtie2samtools       don't cleanup bowtie2samtools output - default false (optional)
      --keep_instrain              don't cleanup instrain output - default false (optional)
      --skip_qc                    skip metawrap qc step - default false (optional)
      --instrain_queue             job queue for instrain - default long (optional)
      --instrain_quick_profile     use quick-profile option for inStrain - default false (optional)
      --bowtie2_samtools_only      only run bowtie2_samtools process - default false (optional)
      -profile                     always use sanger_lsf when running on the farm (mandatory)
      --help                       print this help message (optional)
    """.stripIndent()
}


def VALIDATE_PARAMETERS() {
    // Parameter checking function
    def errors = 0

    if (params.manifest) {
        manifest=file(params.manifest)
        if (!manifest.exists()) {
            log.error("The manifest file specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No manifest file specified. Please specify one using the --manifest option.")
        errors += 1
    }

    if (params.bowtie2_samtools_threads) {
        if (!params.bowtie2_samtools_threads.toString().isInteger()) {
        log.error("Please ensure the bowtie2_samtools_threads parameter is a number")
        errors += 1
        }
    }
    else {
        log.error("Please specify the number of threads for bowtie2samtools using the --bowtie2_samtools_threads option")
        errors += 1
    }

    if (params.instrain_full_output && params.instrain_quick_profile) {
        log.error("the --instrain_full_output and --instrain_quick_profile options are incompatible, please choose one of these options.")
        errors += 1
    }

    if (params.instrain_threads) {
        if (!params.instrain_threads.toString().isInteger()) {
        log.error("Please ensure the instrain_threads parameter is a number")
        errors += 1
        }
    }
    else {
        log.error("Please specify the number of threads for instrain_threads using the --instrain_threads option")
        errors += 1
    }

    if (params.results_dir) {
        results_dir_path=file(params.results_dir)
        if (!results_dir_path.getParent().exists()) {
            log.error("The results directory path specified does not exist.")
            errors += 1
        }
    }
    else {
        log.error("No results directory has been specified, please ensure you provide a value or a default.")
        errors += 1
    }

    if (errors > 0) {
            log.error(String.format("%d errors detected", errors))
            exit 1
        }
}