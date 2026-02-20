#!/usr/bin/env python3


"""
Create a fasta file containing fasta sequences from the sourmash matches
"""

import sys 
import os
import re
import gzip
import shutil

def normalize_ext(ext):
    ext = ext.strip()
    if ext and not ext.startswith("."):
        ext = "." + ext
    return ext

def open_by_ext(path, genomes_file_ext):
    #If they said '.fa.gz' / '.fasta.gz' etc -> use gzip.open
    #Otherwise -> normal open
    if genomes_file_ext.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def subset_fasta(genome_dir, sourmash_genomes, genomes_file_ext):

    # Variable to flag if the directory is not a GCF pattern e.g GCF/001/002/003/

    not_GCF_pattern = True

    # Open the sourmash genomes file and the output subset fasta file
    genome_dir = genome_dir.strip()

    #normalise file extension:
    genomes_file_ext = normalize_ext(genomes_file_ext)

    with open(sourmash_genomes, 'r') as sourmash_file, open('subset_ref_database.fasta', 'ab') as subset_file:
        # Iterate through each genome file name in the sourmash file
        for genome_file_name in sourmash_file: 
            genome_file_name = genome_file_name.strip()

            #blank-line skip
            if not genome_file_name:
                continue
            expected_name = genome_file_name + genomes_file_ext
            # Only execute if directory is in a GCF pattern
            if not_GCF_pattern:
                # Create the path for the genome file based on the GCF pattern
                pattern_re = r"(^[A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3}).*$"
                genome_split_to_path = re.sub(pattern_re, r"database/\1/\2/\3/\4/", genome_file_name)
                genome_file_path = os.path.join(genome_dir, genome_split_to_path, expected_name)
                print("file path: ", genome_file_path)
            if os.path.exists(genome_file_path) and not_GCF_pattern:
                try:
                    with open_by_ext(genome_file_path, genomes_file_ext) as genome_file:
                        shutil.copyfileobj(genome_file, subset_file)
                        continue
                except Exception as e:
                    raise RuntimeError(f"Error processing {genome_file_path}: {e}") from e
            else:
                # If the directory is not in a GCF pattern, search for the genome file in the genome directory
                not_GCF_pattern = False
                print("file path2: ", genome_file_path)

                found = False

                for root, _, files in os.walk(genome_dir):
                    if expected_name in files:
                        genome_file_path = os.path.join(root, expected_name)
                        try:
                            with open_by_ext(genome_file_path, genomes_file_ext) as genome_file:
                                shutil.copyfileobj(genome_file, subset_file)
                                found = True
                                break
                        except Exception as e:
                            raise RuntimeError(f"Error processing {genome_file_path}: {e}") from e

                if not found:
                    raise FileNotFoundError(
                        f"Expected genome file '{expected_name}' not found anywhere under {genome_dir}. "
                        f"Check that your files match --genomes_file_ext exactly."
                    )

if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise TypeError(f"subset_fasta.py takes 3 positional arguments but {len(sys.argv)-1} were given\nUsage: {sys.argv[0]} <genome_dir> <sourmash_genomes> <genomes_file_ext>")
        sys.exit(1)


    genome_dir = sys.argv[1]
    sourmash_genomes = sys.argv[2]
    genomes_file_ext = sys.argv[3]

    subset_fasta(genome_dir, sourmash_genomes, genomes_file_ext)
