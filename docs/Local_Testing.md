# Local Testing <!-- omit in toc -->

- [Pre-requisites](#pre-requisites)
- [`-stub-run` testing](#-stub-run-testing)
  - [Full `-stub-run`](#full--stub-run)
  - [Ascat from counts files `-stub-run`](#ascat-from-counts-files--stub-run)
- [Data based testing](#data-based-testing)
  - [Create a data folder](#create-a-data-folder)
  - [Reference files](#reference-files)
  - [COLO-829 & COLO-829-BL](#colo-829--colo-829-bl)
  - [Full run](#full-run)
  - [Ascat from counts files](#ascat-from-counts-files)

Testing of this workflow can be completed using data provided by the authors of dockstore-cgpwgs.

You can also use `-stub-run` for rapid testing of structural changes with no data requirements.

## Pre-requisites

- Nextflow
- 30GB storage - not needed for `-stub-run`
- samtools - not needed for `-stub-run`

## `-stub-run` testing

Stub run testing takes seconds to complete but only verifies the workflow structure and passing of variables and files.

### Full `-stub-run`

Runs `allele_count` and `ascat` steps:

```bash
rm -rf full_stub_results .nextflow* work
nextflow run main.nf \
    -profile test -stub-run \
    --pairs test_files/colo-cram.csv \
    --outdir full_stub_results
```

### Ascat from counts files `-stub-run`

Only executes `ascat` step:

```bash
rm -rf count_stub_results .nextflow* work
nextflow run main.nf \
    -profile test -stub-run \
    --counts \
    --pairs test_files/colo-counts.csv \
    --outdir count_stub_results
```

## Data based testing

Data based testing takes ~30 minutes once inputs have been generated

### Create a data folder

Create a test folder for the input.  It will need to hold up to ~30GB.

This folder should be used in place of `$TESTDATA` in future commands.

### Reference files

Download and extract the required files:

```
cd $TESTDATA
mkdir GRCh37
curl -sSL ftp.sanger.ac.uk/pub/cancer/dockstore/human/core_ref_GRCh37d5.tar.gz | tar -C GRCh37 --strip-components 1 -zx
curl -sSL ftp.sanger.ac.uk/pub/cancer/dockstore/human/qcGenotype_GRCh37d5.tar.gz | tar -C GRCh37 --strip-components 1 -zx
curl -sSL ftp.sanger.ac.uk/pub/cancer/dockstore/human/CNV_SV_ref_GRCh37d5_brass6+.tar.gz | tar --no-anchor -C GRCh37 --strip-components 1 -zx SnpGcCorrections.tsv
```

### COLO-829 & COLO-829-BL

Low depth copies of this cell line are functional.  Workflow supports CRAM so minimisation of space is considered.

```
cd $TESTDATA

wget http://ngs.sanger.ac.uk/production/cancer/dockstore/cgpwgs/sampled/COLO-829.bam
samtools view --write-index -T GRCh37/genome.fa -Co GRCh37/COLO-829.cram -@ 7 --output-fmt-option normal COLO-829.bam
rm COLO-829.bam

wget http://ngs.sanger.ac.uk/production/cancer/dockstore/cgpwgs/sampled/COLO-829-BL.bam
samtools view --write-index -T GRCh37/genome.fa -Co GRCh37/COLO-829-BL.cram -@ 7 --output-fmt-option normal COLO-829-BL.bam
rm COLO-829-BL.bam
```

### Full run

Executes the `allele_count` and `ascat` steps.

1. Copy the csv file `cp test_data/colo-cram.csv $TESTDATA/.`
1. Update `colo-cram.csv` replacing `TESTDATA` with the expanded value (4 instances)
1. Run following commands:

```bash
rm -rf full_results .nextflow* work
nextflow run main.nf \
    -profile test \
    --cpus_counts 3 \
    --genome_fa $TESTDATA/GRCh37/genome.fa \
    --snp_gc_corr $TESTDATA/GRCh37/ascat/SnpGcCorrections.tsv \
    --gender_snps $TESTDATA/GRCh37/gender.tsv \
    --pairs $TESTDATA/colo-cram.csv \
    --outdir $TESTDATA/full_results
```

### Ascat from counts files

Outputs from "Full run" required

1. Copy the csv file `cp test_data/colo-counts.csv $TESTDATA/.`
1. Update `colo-counts.csv` replacing `TESTDATA` with the expanded value (6 instances)
1. Run following commands:

```bash
rm -rf count_results .nextflow* work
nextflow run main.nf \
    -profile test \
    --counts \
    --genome_fa $TESTDATA/GRCh37/genome.fa \
    --snp_gc_corr $TESTDATA/GRCh37/ascat/SnpGcCorrections.tsv \
    --gender_snps $TESTDATA/GRCh37/gender.tsv \
    --pairs $TESTDATA/colo-cram.csv \
    --outdir $TESTDATA/count_results
```
