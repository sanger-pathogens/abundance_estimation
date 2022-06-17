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
  --skip_qc                    skip metawrap qc step - default false, use true if qc is not needed
  --cleanup                    cleanup intermediate files - default true, use false if cleanup is not wanted
  -profile                     always use sanger_lsf when running on the farm
```

## Generating manifests
### Step 1:
Obtain fastq paths:
`ls -d -1 <path>/*.fastq.gz > fastq_paths.txt`
### Step 2:
Generate manifest:
`./generate_manifest.sh fastq_paths.txt`

This will output the manifest to `manifest.csv` which can be fed into the nextflow pipeline