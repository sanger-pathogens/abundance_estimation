#!/usr/bin/env bash

path_to_taxonomy_lookup="/data/pam/team230/wr7/scratch/rvi/instrain/dbs/viral_refseq/2023_09_09/viral.1.1.genomic_custom.stb"

declare -A seq_id_to_name; while read key value; do seq_id_to_name["$key"]="$value"; done < "${path_to_taxonomy_lookup}"

for file in *.tsv; do
  paste -d '\t' <(echo genome; while read key; do echo ${seq_id_to_name["$key"]};done < <(tail +2 "$file" | sed "s/.*viral.1.1.genomic_sorted.id_\(.*\).fna.*/\1/g")) <(awk -F $'\t' '{for (i=2; i<NF; i++) printf $i "\t"; print $NF}' "$file") > "${file%.tsv}_fixed.tsv"
done
