# Amazon AWS cloud

You can run this pipeline on AWS cloud.

You need to set up AWS batch first, as decribed [here](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/tipi).

Then, copy `nextflow.aws.config.template` to `nextflow.aws.config` and modify the `aws` values to match your setup (if you do not want to put access keys into your config file, you can also use environment variables).

Finally, upload your data to AWS. You can use S3 (cheaper) or EFS for `stb_file`, `genome_dir`, `outdir`, and `workdir`; for `sourmash_db` and `bmtagger_db` you need EFS. `manifest` can stay on your local disk.

Your input files (eg `.bam`) need to be copied into the cloud as well, best to S3. You will have to change your manifest files accordingly, see [example](manifest.aws.csv.template).

You can now run your pipeline from your local machine, like so:

```bash
nextflow -C nextflow.config.aws run main.nf --manifest mainfest.csv
```
