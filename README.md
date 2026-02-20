# Abundance estimation nextflow pipeline

## Usage

```bash
nextflow run main.nf \
  --manifest <manifest.csv> \
  --outdir ./results \
  --instrain_full_output false \
  --cleanup_intermediate_files false \
  --skip_qc false \
  --stb_file /data/pam/software/GTDB/gtdb_genomes_reps_r226.stb \
  --genome_dir /data/pam/software/GTDB/release226/genomic_files_reps/gtdb_genomes_reps_r226 \
  --sourmash_db /data/pam/software/sourmash/GTDB/release226/gtdb-rs226-reps.k31.sig.zip \
  --genomes_file_ext .fasta \
  --instrain_quick_profile false \
  --bowtie2_samtools_only false
```

Note: parameters marked with (optional*) are only optional if you use the default GTDB resources.
If you supply a custom `--sourmash_db`, you must also provide matching `--genome_dir`, `--genomes_file_ext`,
and `--stb_file` built from the same reference set.

### Parameters

- `--manifest` Manifest containing paths to fastq files with headers `ID,R1,R2`. (mandatory)
- `--outdir` Name of results folder. [default: `./results`] (optional)
- `--instrain_full_output` Get full inStrain output. [default: false] (optional)
- `--cleanup_intermediate_files` Cleanup intermediate files. [default: false] (optional)
- `--skip_qc` Skip metawrap qc step. [default: false] (optional)
- `--stb_file` Supply stb file. [default: `/data/pam/software/GTDB/gtdb_genomes_reps_r226.stb`] (optional*)
- `--genome_dir` Supply genome folder. [default: `/data/pam/software/GTDB/release226/genomic_files_reps/gtdb_genomes_reps_r226`] (optional*)
- `--sourmash_db` Supply sourmash database. [default: `/data/pam/software/sourmash/GTDB/release226/gtdb-rs226-reps.k31.sig.zip`] (optional*)
- `--genomes_file_ext` File extension for reference genomes (e.g. `.fna`, `.fasta`, `_genomic.fna.gz`). [default: `.fasta`] (optional*)
- `--instrain_quick_profile` Use quick-profile option for inStrain. [default: false] (optional)
- `--bowtie2_samtools_only` Only run bowtie2_samtools process. [default: false] (optional)
- `--help` Print this help message. (optional)

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

```bash
nextflow/22.10
ISG/singularity/3.6.4
```

## Resource requirements

`bowtie2samtools` - 250GB ( ~2hr) submits with 250Gb and then on retry escalate to 350Gb, 8 CPUs

`inStrain` - 300 GB ( ~ 12hr) submits with 300Gb then on retry escalates to 400Gb, 8 CPUs

It’s not particularly uncommon for these things to run out of memory, so they do need to retry on a few samples for each run

## Custom sourmash database generation

To begin you need to have a set of input genomic sequences:

```
GCA_900538355.1.fasta
GCA_900549805.1.fasta
GCA_900550455.1.fasta
```

In this example the extension is `.fasta`. The same steps apply to `.fna` files.

To generate a sourmash database from these files, you first need to produce a sketch for each input fasta:

```bash
module load sourmash/4.5.0--hdfd78af_0
sourmash sketch dna -p scaled=1000,k=31 db/*.fasta
```

This will produce a collection of signature files ending in `.sig`:

```
GCA_900538355.1.fasta.sig
GCA_900549805.1.fasta.sig
GCA_900550455.1.fasta.sig
```

This pipeline requires that the name of the signature within the signature file is the same as the basename of the file i.e.

```
GCA_900538355.1
```

This can be produced using the command `sourmash signature rename` included in the sourmash package.

To rename a collection of signature file you can use the following commands. First list the files into a list to use to rename:

```bash
ls *.sig > filelist
```

Then run a loop over this list using the sourmash script to rename the signature to the filename:

```bash
cat filelist | while read line;
do
    NAME=$(basename "$line" .fasta.sig)
    echo $NAME
    sourmash signature rename $NAME.fasta.sig "$NAME" -o $NAME.sig
done
```

Cleanup the original signatures if you only want the renamed files:

```bash
rm -f $(cat filelist)
```

Once complete, index the output signature files into an indexed record:

```bash
sourmash index -k 31 all-genomes *.sig
```

Supply the index as an argument to the pipeline option `--sourmash_db`.

In this scenario you will also need to include the following option (specifying the file extension of the input sequences that were used to build the sourmash index):

```bash
--genomes_file_ext .fasta
```

And point to the genome dir where the .fasta files are stored:

```bash
--genome_dir <path_to_fasta_files>
```

You must also provide an STB file built from the same reference genomes.

Build an STB mapping file from reference FASTA headers:

```bash
cat > /tmp/make_stb.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

for f in /path/to/refs/*.fasta; do
  bin=$(basename "$f" .fasta)
  awk -v b="$bin" '
    /^>/ {
      header=$0
      sub(/^>/,"",header)
      split(header, a, /[ \t]/)
      print a[1] "\t" b
    }
  ' "$f"
done
EOF

bash /tmp/make_stb.sh > /path/to/refs/custom.stb
```

Then pass it to the pipeline:

```bash
--stb_file /path/to/refs/custom.stb
```

Make sure the `scaled` parameter used to build the sourmash database matches the pipeline sketch settings.
If they do not match, sourmash gather may return few or no hits.


## Amazon AWS

See [the AWS README](README.aws.md).
