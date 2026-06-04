# abundance_estimation

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[[_TOC_]]

## Pipeline overview

**abundance_estimation** is a Nextflow DSL2 pipeline for estimating the relative abundance of microbial species from metagenomic short-read data. It uses [Sourmash](https://sourmash.readthedocs.io/) for rapid genome identification, [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) for competitive read mapping to reference genomes, and [inStrain](https://instrain.readthedocs.io/) for strain-resolved abundance profiling.

The pipeline performs the following steps:

1. **QC** (optional) — adapter trimming and host read removal with MetaWrap QC (TrimGalore + BMTagger).
2. **Read merging** — per-sample FASTQ lanes are merged.
3. **Reference selection** — Sourmash sketches the reads and queries the GTDB genome database to identify the subset of reference genomes present in each sample.
4. **Competitive mapping** — Bowtie2 indexes the selected reference genomes and maps reads competitively to them.
5. **Abundance profiling** — inStrain profiles strain-level abundance from the mapping output.

## Usage

### Quickstart

#### From source code

1. Clone this repository (including submodules):

   ```bash
   git clone --recurse-submodules https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/abundance_estimation.git
   cd abundance_estimation
   ```

2. Run with `singularity`:

   ```bash
   nextflow run main.nf \
       -profile singularity \
       --manifest manifest.csv \
       --outdir my_output
   ```

3. Once the run has finished, clean up intermediate files:

   ```bash
   rm -rf work .nextflow*
   ```

#### Using on the Sanger farm

Load Nextflow and Singularity:

```bash
module load nextflow ISG/singularity
```

Submit to LSF:

```bash
bsub -o output.o -e error.e -q oversubscribed -R "select[mem>4000] rusage[mem=4000]" -M4000 \
    nextflow run main.nf \
        --manifest manifest.csv \
        --outdir my_output
```

### Input

#### Manifest (`--manifest`)

A CSV file with the required header `ID,R1,R2`, containing per-sample paths to paired `.fastq.gz` files:

```
ID,R1,R2
sampleA,/path/to/sampleA_1.fastq.gz,/path/to/sampleA_2.fastq.gz
sampleB,/path/to/sampleB_1.fastq.gz,/path/to/sampleB_2.fastq.gz
```

#### Generating a manifest

**Sanger users:** the [manifest_generator](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/manifest_generator/) tool can generate a compatible `ID,R1,R2` manifest from a directory of FASTQ files or from iRODS.

### Output

Results are written to `--outdir` (default: `./results`):

```
results/
  <sample_ID>/
    instrain/                        # inStrain profiling output
      output/
        <sample_ID>_genome_info.tsv  # Per-genome abundance and coverage
        <sample_ID>_mapping_info.tsv # Per-read mapping details
  bowtie2/
    <sample_ID>.bam                  # Sorted BAM file (competitive mapping)
    <sample_ID>.overall_mapping_rate.txt
  sourmash/
    <sample_ID>_gather.csv           # Sourmash gather results (genome matches)
```

### Parameters

| Option                         | Type      | Default                                                                        | Description                                                        |
| ------------------------------ | --------- | ------------------------------------------------------------------------------ | ------------------------------------------------------------------ |
| `--manifest`                   | `path`    | (required)                                                                     | Input manifest CSV with header `ID,R1,R2`.                         |
| `--outdir`                     | `path`    | `./results`                                                                    | Directory where results are written.                               |
| `--skip_qc`                    | `boolean` | `false`                                                                        | Skip MetaWrap QC (adapter trimming and host read removal).         |
| `--stb_file`                   | `path`    | `/data/pam/software/GTDB/gtdb_genomes_reps_r226.stb`                           | Sample-to-bin (STB) mapping file for inStrain.                     |
| `--genome_dir`                 | `path`    | `/data/pam/software/GTDB/release226/genomic_files_reps/gtdb_genomes_reps_r226` | Directory containing GTDB reference genome FASTAs.                 |
| `--sourmash_db`                | `path`    | `/data/pam/software/sourmash/signatures_zipped/gtdb_genomes_reps_r220.zip`     | Sourmash genome signature database.                                |
| `--instrain_full_output`       | `boolean` | `false`                                                                        | Publish full inStrain output (large).                              |
| `--instrain_quick_profile`     | `boolean` | `false`                                                                        | Use inStrain `quick_profile` mode (faster, less detail).           |
| `--bowtie2_samtools_only`      | `boolean` | `false`                                                                        | Run only Bowtie2 mapping and Samtools steps, skipping inStrain.    |
| `--cleanup_intermediate_files` | `boolean` | `false`                                                                        | Delete intermediate files (trimmed FASTQs, sorted BAMs) after use. |

### Advanced usage

#### Mapping only (no inStrain)

To generate BAM files without running inStrain (useful for downstream analysis or to reduce compute):

```bash
nextflow run main.nf --manifest manifest.csv --bowtie2_samtools_only true --outdir my_output
```

#### Skipping host removal and trimming

If reads have already been quality-controlled:

```bash
nextflow run main.nf --manifest manifest.csv --skip_qc true --outdir my_output
```

#### Using a custom Sourmash database

Supply your own Sourmash signature database and matching genome directory:

```bash
nextflow run main.nf \
    --manifest manifest.csv \
    --sourmash_db /path/to/custom.sig.zip \
    --genome_dir /path/to/genomes \
    --stb_file /path/to/custom.stb \
    --outdir my_output
```

See the existing README for instructions on building a custom Sourmash database.

### Dependencies

All software dependencies are containerised. The following databases must be available (Sanger HPC defaults are pre-configured):

- **Sourmash database** (`--sourmash_db`): GTDB genome signature zip file.
- **Reference genome directory** (`--genome_dir`): directory of GTDB reference FASTA files.
- **STB file** (`--stb_file`): sample-to-bin mapping file required by inStrain.
- **BMTagger database**: required when `--skip_qc false` (host read removal). Pre-configured on the Sanger HPC.

**Resource requirements**: Bowtie2 mapping requires ~250 GB RAM; inStrain requires ~300 GB RAM. These processes will automatically retry with increased memory on failure.

## Software versions

| Software           | Version | Image                                              |
| ------------------ | ------- | -------------------------------------------------- |
| Sourmash           | 4.5.0   | `quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0` |
| Bowtie2 + Samtools | —       | `quay.io/sangerpathogens/bowtie2-samtools:1.1-c1`  |
| inStrain           | 1.9.0   | `quay.io/sangerpathogens/instrain:1.9.0`           |

See `modules/` for pinned container versions.

## Troubleshooting

- **Out-of-memory errors**: Bowtie2 and inStrain are very memory-intensive. The pipeline is configured to retry with more memory on failure. Ensure the HPC queue has nodes with sufficient RAM (>300 GB).
- **Sourmash finds no matches**: ensure `--sourmash_db` and `--genome_dir` are consistent (same GTDB release). Check that reads are of sufficient quality and depth.
- **Resuming a failed run**: add `-resume` to restart from cached intermediate results. Note: if `--cleanup_intermediate_files true` is set, files deleted in earlier runs cannot be reused.
- For further help, check `.nextflow.log` and the per-process logs in the `work/` directory.

## Issues and Contributions

**GitHub users:** if you find an issue with this pipeline, or would like to suggest an improvement, please log an issue or open a pull request on this repository.

**Sanger users:** if you need internal support, you can raise an issue on the PAM Freshservice portal: https://sanger.freshservice.com/support/catalog/items/426
