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
    genome_dir = genome_dir.strip()
    with open(sourmash_genomes, 'r') as sourmash_file, open('subset_ref_database.fasta', 'ab') as subset_file:
        # Iterate through each genome file name in the sourmash file
        for genome_file_name in sourmash_file: 
            genome_file_name = genome_file_name.strip()
            # Only execute if directory is in a GCF pattern
            if not_GCF_pattern:
                # Create the path for the genome file based on the GCF pattern
                pattern_re = r"(^[A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3}).*$"
                genome_split_to_path = re.sub(pattern_re, r"database/\1/\2/\3/\4/", genome_file_name)
                genome_file_path = os.path.join(genome_dir, genome_split_to_path, genome_file_name)
                print("file path: ", genome_file_path)
            if os.path.exists(genome_file_path) and not_GCF_pattern:
                try:
                    with gzip.open(genome_file_path, 'rb') as genome_file:
                        shutil.copyfileobj(genome_file, subset_file)
                except FileNotFoundError as e:
                    print(f"Error {genome_file_path} not found: {e}")
                except BaseException as e:
                    print(f"Error processing {genome_file_path}: {e}")
            else:
                # If the directory is not in a GCF pattern, search for the genome file in the genome directory
                not_GCF_pattern = False
                print("file path2: ", genome_file_path)
                for root, _, files in os.walk(genome_dir):
                    if genome_file_name in files:
                        genome_file_path = os.path.join(root, genome_file_name)
                        try:
                            with gzip.open(genome_file_path, 'rb') as genome_file:
                                shutil.copyfileobj(genome_file, subset_file)
                        except FileNotFoundError as e:
                            print(f"Error {genome_file_path} not found: {e}")
                        except BaseException as e:
                            print(f"Error processing {genome_file_path}: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise TypeError(f"two() takes 2 positional arguments but {len(sys.argv)-1} were given\nUsage: {sys.argv[0]} <genome_dir> <sourmash_genomes>")
        sys.exit(1)


    genome_dir = sys.argv[1]
    sourmash_genomes = sys.argv[2]

    subset_fasta(genome_dir, sourmash_genomes)