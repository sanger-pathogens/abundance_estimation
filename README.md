# Abundance estimation nextflow pipeline

## Usage
```
nextflow run .
  --manifest                   Manifest containing paths to fastq files (mandatory)
  --btidx                      bowtie index - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.bt2 (optional)               
  --genome_file                genome file - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.fasta (optional)
  --stb_file                   stb file - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.stb (optional)
  --bowtie2_samtools_threads   threads - default: 8 (optional)
  --instrain_threads           threads - default: 8 (optional)
  --instrain_full_output       get full instrain output - default false (optional)
  --keep_metawrap_qc           don't cleanup metawrap_qc output - default false (optional)
  --keep_bowtie2samtools       don't cleanup bowtie2samtools output - default false (optional)
  --keep_instrain              don't cleanup instrain output - default false (optional)
  --skip_qc                    skip metawrap qc step - default false (optional)
  --instrain_queue             job queue for instrain - default long (optional)
  --instrain_quick_profile     use quick-profile option for inStrain - default false (optional)
  --bowtie2_samtools_only      only run bowtie2_samtools process - default false (optional)
  -profile                     always use sanger_lsf when running on the farm (mandatory)
  --help                       print this help message (optional)
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
