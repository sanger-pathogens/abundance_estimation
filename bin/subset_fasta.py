#!/usr/bin/env python3


"""
Create a fasta file containing fasta sequences from the sourmash matches
"""

import sys 
import os
import re
import gzip
import shutil

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def normalize_ext(ext):
    ext = ext.strip()
    if ext and not ext.startswith("."):
        ext = "." + ext
    return ext


def filename_variants(name, ext):
    if ext and not name.endswith(ext):
        return [name + ext, name]
    return [name]


def subset_fasta(genome_dir, sourmash_genomes, genomes_file_ext):

    # Variable to flag if the directory is not a GCF pattern e.g GCF/001/002/003/

    not_GCF_pattern = True

    # Open the sourmash genomes file and the output subset fasta file
    genome_dir = genome_dir.strip()
    genomes_file_ext = normalize_ext(genomes_file_ext)
    with open(sourmash_genomes, 'r') as sourmash_file, open('subset_ref_database.fasta', 'ab') as subset_file:
        # Iterate through each genome file name in the sourmash file
        for genome_file_name in sourmash_file: 
            genome_file_name = genome_file_name.strip()
            genome_file_path = None

            # Only execute if directory is in a GCF pattern
            if not_GCF_pattern:
                # Create the path for the genome file based on the GCF pattern
                pattern_re = r"(^[A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3}).*$"
                if re.match(pattern_re, genome_file_name):
                    genome_split_to_path = re.sub(pattern_re, r"database/\1/\2/\3/\4/", genome_file_name)
                    for name in filename_variants(genome_file_name, genomes_file_ext):
                        genome_file_path = os.path.join(genome_dir, genome_split_to_path, name)
                        print("file path: ", genome_file_path)
                        if os.path.exists(genome_file_path):
                            try:
                                with open_maybe_gzip(genome_file_path) as genome_file:
                                    shutil.copyfileobj(genome_file, subset_file)
                            except FileNotFoundError as e:
                                print(f"Error {genome_file_path} not found: {e}")
                            except BaseException as e:
                                print(f"Error processing {genome_file_path}: {e}")
                            continue
                        if os.path.exists(genome_file_path + ".gz"):
                            try:
                                with open_maybe_gzip(genome_file_path + ".gz") as genome_file:
                                    shutil.copyfileobj(genome_file, subset_file)
                            except FileNotFoundError as e:
                                print(f"Error {genome_file_path}.gz not found: {e}")
                            except BaseException as e:
                                print(f"Error processing {genome_file_path}.gz: {e}")
                            continue
                    else:
                        not_GCF_pattern = False
                else:
                    not_GCF_pattern = False

            if not_GCF_pattern:
                continue
            else:
                # If the directory is not in a GCF pattern, search for the genome file in the genome directory
                for root, _, files in os.walk(genome_dir):
                    for name in filename_variants(genome_file_name, genomes_file_ext):
                        if name in files:
                            genome_file_path = os.path.join(root, name)
                            try:
                                with open_maybe_gzip(genome_file_path) as genome_file:
                                    shutil.copyfileobj(genome_file, subset_file)
                            except FileNotFoundError as e:
                                print(f"Error {genome_file_path} not found: {e}")
                            except BaseException as e:
                                print(f"Error processing {genome_file_path}: {e}")
                            break
                        if name + ".gz" in files:
                            genome_file_path = os.path.join(root, name + ".gz")
                            try:
                                with open_maybe_gzip(genome_file_path) as genome_file:
                                    shutil.copyfileobj(genome_file, subset_file)
                            except FileNotFoundError as e:
                                print(f"Error {genome_file_path} not found: {e}")
                            except BaseException as e:
                                print(f"Error processing {genome_file_path}: {e}")
                            break



if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise TypeError(f"two() takes 3 positional arguments but {len(sys.argv)-1} were given\nUsage: {sys.argv[0]} <genome_dir> <sourmash_genomes> <genomes_file_ext>")
        sys.exit(1)


    genome_dir = sys.argv[1]
    sourmash_genomes = sys.argv[2]
    genomes_file_ext = sys.argv[3]

    subset_fasta(genome_dir, sourmash_genomes, genomes_file_ext)
