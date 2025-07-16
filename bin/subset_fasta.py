#!/usr/bin/env python3


"""
Create a fasta file containing fasta sequences from the sourmash matches
"""

import sys 
import os
import re
import gzip
import shutil


def subset_fasta(genome_dir, sourmash_genomes):

    # Variable to flag if the directory is not a GCF pattern e.g GCF/001/002/003/

    not_GCF_pattern = True

    # Open the sourmash genomes file and the output subset fasta file
    with open(sourmash_genomes, 'r') as sourmash_file:
        with open('subset_ref_database.fasta', 'ab') as subset_file:
            # Iterate through each genome file name in the sourmash file
            for genome_file_name in sourmash_file: 
                # Only execute if directory is in a GCF pattern
                if not_GCF_pattern:
                    # Create the path for the genome file based on the GCF pattern
                    genome_split_to_path = re.sub(r"(^[A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3}).*$",r"database/\1/\2/\3/\4/", genome_file_name)
                    genome_file_path = os.path.join(genome_dir.strip(), genome_split_to_path.strip(), genome_file_name.strip())
                    print("file path: ", genome_file_path)
                if os.path.exists(genome_file_path) and not_GCF_pattern:
                    with gzip.open (genome_file_path, 'rb') as genome_file:
                        shutil.copyfileobj(genome_file, subset_file)
                else:
                    # If the directory is not in a GCF pattern, search for the genome file in the genome directory
                    not_GCF_pattern = False
                    print("file path2: ", genome_file_path)
                    for root, _, files in os.walk(genome_dir.strip()):
                        if genome_file_name.strip() in files:
                            genome_file_path = os.path.join(root, genome_file_name.strip())
                            with gzip.open(genome_file_path, 'rb') as genome_file:
                                shutil.copyfileobj(genome_file, subset_file)
                            break

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python subset_fasta.py <genome_dir> <sourmash_genomes>")
        sys.exit(1)

    genome_dir = sys.argv[1]
    sourmash_genomes = sys.argv[2]

    subset_fasta(genome_dir, sourmash_genomes)