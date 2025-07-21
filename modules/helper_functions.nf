def validate_path_param(
    param_option,
    param,
    type="file",
    mandatory=true) {
        if (param.startsWith("/efs/")) { // This is a hack to avoid checking the existence of files in the EFS
            return 0
        }
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

def validate_outdir(outdir) {
    outdir = file(outdir)
    if (outdir.exists() && !outdir.isDirectory()) {
        log.error("The given outdir '${outdir}' is not a directory.")
        return 1
    }
    return 0
}

def validate_parameters() {
    def errors = 0

    errors += validate_path_param("--manifest", params.manifest)
    errors += validate_path_param("--stb_file", params.stb_file)
    errors += validate_path_param("--sourmash_db", params.sourmash_db)
    errors += validate_path_param("--genome_dir", params.genome_dir, type="directory")
    errors += validate_number_param("--queue_size", params.queue_size)
    errors += validate_outdir(params.outdir)

    if (params.instrain_full_output && params.instrain_quick_profile) {
        log.error("the --instrain_full_output and --instrain_quick_profile options are incompatible, please choose one of these options.")
        errors += 1
    }

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}