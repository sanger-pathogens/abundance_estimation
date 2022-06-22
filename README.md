# Abundance estimation nextflow pipeline

## Usage
```
nextflow run abundance_estimation.nf
  --manifest                   Manifest containing paths to fastq files (mandatory)
  --btidx                      bowtie index - default: /lustre/scratch118/infgen/team162/shared/gtdb_genomes_reps_r202/gtdb_genomes_reps_r202.fasta.bt2 (optional)               
  --genome_file                genome file - default: /lustre/scratch118/infgen/team162/shared/gtdb_genomes_reps_r202/gtdb_genomes_reps_r202.fasta (optional)
  --stb_file                   stb file - default: /lustre/scratch118/infgen/team162/shared/gtdb_genomes_reps_r202/gtdb_genomes_reps_r202.stb
  --bowtie2_samtools_threads   threads - default: 16 (optional)
  --instrain_threads           threads - default: 16 (optional)
  --full_output                get full instrain output - default false, use true if full output required (optional)
  --skip_qc                    skip metawrap qc step - default false, use true if qc is not needed (optional)
  --no_cleanup                 don't cleanup intermediate files - default false (optional)
  -profile                     always use sanger_lsf when running on the farm (mandatory)
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
