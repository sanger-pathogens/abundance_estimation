#!/bin/bash
sourmash_genomes="${1}"
genome_dir="${2}"
genomes_file_ext="${3}"
while read genome
    do
      genome_split_to_path=$(echo ${genome} | sed -e 's|\(GC[AF]\)_\([0-9]\{3\}\)\([0-9]\{3\}\)\([0-9]\{3\}\)\.[0-9]|\1/\2/\3/\4|')
      genome_file_path="${genome_dir}/${genome_split_to_path}/${genome}${genomes_file_ext}"
      [ -e ${genome_file_path} ] || genome_file_path="${genome_dir}/${genome}${genomes_file_ext}"
      zcat -f ${genome_file_path} >> subset_ref_database.fasta
done < ${sourmash_genomes}