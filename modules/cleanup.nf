process cleanup_sorted_bam_files {
    /**
    * Cleanup intermediate files
    */

    input:
        file(sorted_bam_file)

    script:
        """
        # Remove sorted bam files
        orig_bam=\$(readlink -f ${sorted_bam_file})
        rm ${sorted_bam_file} # Remove symlink
        rm \${orig_bam} # Remove original files
        """
}

process cleanup_trimmed_fastq_files {
    /**
    * Cleanup intermediate files
    */

    input:
        tuple val(sample_id), file(first_read), file(second_read)

    script:
        """
        # Remove trimmed fastq files
        orig_fastqs=\$(readlink -f ${first_read} ${second_read})
        rm ${first_read} ${second_read} # Remove symlink
        rm \${orig_fastqs} # Remove original files
        """
}

process cleanup_instrain_output {
    /**
    * Cleanup unused output
    */

    input:
         path(workdir)
         val(sample_id)
    script:
        """
        # Remove instrain results
        instrain_dir=\$(cat $workdir)
        cd \$instrain_dir
        rm -rf $sample_id*
        """
}
