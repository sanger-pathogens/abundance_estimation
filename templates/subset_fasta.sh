#!/bin/bash
# shellcheck disable=SC1083,SC2034,SC2154
# above directive prevents shellcheck to report use of
# SC1083: litteral syntax \${}, required as this is a template script
# SC2034: apparently unused variables, which derives from the above
# SC2154: variables referenced but apparently not assigned; this results from nextflow variables being injected into the bash code
while read genome
    do
      genome_split_to_path=`echo "\${genome}" | sed -e 's|\\(GC[AF]\\)_\\([0-9]\\{3\\}\\)\\([0-9]\\{3\\}\\)\\([0-9]\\{3\\}\\)\\.[0-9]|\\1/\\2/\\3/\\4|'`
      genome_file_path="${genome_dir}/\${genome_split_to_path}/\${genome}${genomes_file_ext}"
      [ -e \${genome_file_path} ] || genome_file_path="${genome_dir}/\${genome}${genomes_file_ext}"

      if [[ \${genome_file_path} == s3:* ]] # S3 mount
      then
        ${aws_cli} s3 --no-progress cp \${genome_file_path} - | zcat >> subset_ref_database.fasta
      else # File system
        zcat -f \${genome_file_path} >> subset_ref_database.fasta
      fi

done < ${sourmash_genomes}