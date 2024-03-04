# Abundance estimation nextflow pipeline

## Usage
```
nextflow run main.nf
      --manifest                      Manifest containing paths to fastq files. with headers ID,R1,R2. (mandatory)
      --outdir                        Name of results folder. [default: ./results] (optional)
      --bowtie2_samtools_threads      Threads for bowtie2 and samtools. [default: 4] (optional)
      --instrain_threads              Threads for instrain. [default: 4] (optional)
      --instrain_full_output          Get full instrain output. [default: false] (optional)
      --cleanup_intermediate_files    Cleanup intermediate files. [default: false] (optional)
      --skip_qc                       Skip metawrap qc step. [default: false] (optional)
      --stb_file                      Supply stb file. [default: /lustre/scratch125/pam/data/software/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.stb] (optional)
      --genome_dir                    Supply genome folder. [default: /data/pam/team162/shared/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207_genome_dir] (optional)
      --sourmash_db                   Supply sourmash database. [default: /data/pam/team162/shared/sourmash_db/gtdb-rs207.genomic-reps.dna.k31.zip] (optional)
      --instrain_quick_profile        Use quick-profile option for inStrain. [default: false] (optional)
      --bowtie2_samtools_only         Only run bowtie2_samtools process. [default: false] (optional)
      --help                          Print this help message. (optional)
```

## Generating manifests

If your data is stored in the PaM informatics pipeline system, you can use the following method:

`./generate_manifest_from_lanes.sh -l <lanes_file>`

For more information, run:
`./generate_manifest_from_lanes.sh -h`

If your data is not stored in the PaM informatics pipeline system, use the following method:
### Step 1:
Obtain fastq paths:
`ls -d -1 <path>/*.fastq.gz > fastq_paths.txt`
### Step 2:
Generate manifest:
`./generate_manifest.sh fastq_paths.txt`

This will output the manifest to `manifest.csv` which can be fed into the nextflow pipeline

## Development
For development, smaller test databases are available, these will significantly reduce the run time and resource requirements:
8 CPUs and 50GB memory will be sufficient

## Dependencies
This pipeline relies on the following modules:
```
nextflow/22.10
ISG/singularity/3.6.4
```

## Resource requirements
`bowtie2samtools` - 250GB ( ~2hr) submits with 250Gb and then on retry escalate to 350Gb, 8 CPUs

`inStrain` - 300 GB ( ~ 12hr) submits with 300Gb then on retry escalates to 400Gb, 8 CPUs

It’s not particularly uncommon for these things to run out of memory, so they do need to retry on a few samples for each run
