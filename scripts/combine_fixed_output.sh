#!/usr/bin/env bash

set -eo pipefail

readarray -t species_array < <(tail -q -n +2 *_fixed.tsv | awk -F '\t' '{print $1}' | sort -u)

for file in *_fixed.tsv; do
  unset lookup
  declare -A lookup
  while IFS=$'\t' read species filtered_reads; do
    lookup["$species"]=$filtered_reads
  done < <(tail +2 $file | awk -F $'\t' -v OFS=$'\t' '{print $1,$(NF-4)}')
  # For every elem in the species array, lookup the correponding number of filtered reads (empty if no species present)
  echo ${file%_genome_info_fixed.tsv} > ${file}_species_lookup.tmp
  for species in "${species_array[@]}"; do
    reads=${lookup["$species"]}
    if [ -z "$reads" ]; then reads=0; fi
    echo $reads
  done >> ${file}_species_lookup.tmp
done

cat <(echo) <(printf '%s\n' "${species_array[@]}") > row_labels.tmp
paste row_labels.tmp *_species_lookup.tmp > summary.tsv
rm *.tmp 
