# Abundance estimation nextflow pipeline

## Usage
```
nextflow run abundance_estimation.nf
  --manifest                   Manifest containing paths to fastq files (mandatory)
  --btidx                      bowtie index - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.bt2 (optional)               
  --genome_file                genome file - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.fasta (optional)
  --stb_file                   stb file - default: /lustre/scratch125/pam/pathogen/pathpipe/gtdb/gtdb_genomes_reps_r207/gtdb_genomes_reps_r207.stb
  --bowtie2_samtools_threads   threads - default: 16 (optional)
  --instrain_threads           threads - default: 16 (optional)
  --instrain_full_output       get full instrain output - default false (optional)
  --keep_metawrap_qc           don't cleanup metawrap_qc output - default false (optional)
  --keep_bowtie2samtools       don't cleanup bowtie2samtools output - default false (optional)
  --keep_instrain              don't cleanup instrain output - default false (optional)
  --skip_qc                    skip metawrap qc step - default false (optional)
  --instrain_queue             job queue for instrain - default normal (optional)
  --bowtie2samtools_queue      job queue for bowtie2samtools - default normal (optional)
  -profile                     always use sanger_lsf when running on the farm (mandatory)
  --help                       print this help message (optional)
```

## Generating manifests
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
bowtie2/2.3.5--py37he860b03_0
instrain/1.5.4
samtools/1.9
metawrap_custom/1.3.2-c11
nextflow/22.04.5-5708
```
