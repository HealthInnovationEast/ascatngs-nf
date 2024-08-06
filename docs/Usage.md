# Usage <!-- omit in toc -->

- [Overview](#overview)
- [Required args](#required-args)
  - [`--pairs`](#--pairs)
  - [`--genome_fa`](#--genome_fa)
  - [`--snp_gc_corr`](#--snp_gc_corr)
  - [`--gender_snps`](#--gender_snps)
  - [`--counts`](#--counts)
- [Optional args](#optional-args)
  - [Max job request options](#max-job-request-options)

## Overview

The workflow can be used in 2 ways:

1. BAM/CRAM sample inputs

- Generates allele counts and sex verification for each sample (parallel processes)
- Runs ASCAT

2. `count.gz` + `is_male.txt` sample inputs

- Runs ASCAT only

Option (2) is intended for refinement of solution by providing refined purity & ploidy values.

Both routes result in the same set of output files.

## Required args

Workflow would fail if these are not defined, further details below table:

| Parameter     | Description                                                                                                              | Type      | Default              | Required | Hidden |
| ------------- | ------------------------------------------------------------------------------------------------------------------------ | --------- | -------------------- | -------- | ------ |
| `pairs`       | CSV file of input info, see docs/Usage                                                                                   | `string`  | pairs.csv            | True     |        |
| `genome_fa`   | Genome fasta (expects co-located index)                                                                                  | `string`  | genome.fa            | True     |        |
| `snp_gc_corr` | SNP loci and GC corrections file                                                                                         | `string`  | SnpGcCorrections.tsv | True     |        |
| `gender_snps` | SNPs specific to Y chromosome                                                                                            | `string`  | gender.tsv           | True     |        |
| `counts`      |                                                                                                                          | `boolean` |                      |          |        |
| `cpus_counts` | Increasing can reduce run time, however values over 4 are not recommended                                                | `integer` | 4                    |          |        |
| `outdir`      | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string`  | results              |          |        |

### `--pairs`

CSV file with header of:

```
groupId,sampleId,type,protocol,platform,reads,readIdx,purity,ploidy,isMale
```

- groupId = Links tumour & normal samples
- sampleId = Name of the sample
- type = case/control
- protocol = e.g. WGS, WES
- platform = e.g. ILLUMINA
- reads = BAM/CRAM file or count.gz (see --counts option)
- readIdx = BAM/CRAM index file or count.gz (see --counts option)
- purity = NA or Purity (rho) setting for manual setting of sunrise plot location
- ploidy = NA or Ploidy (psi) setting for manual setting of sunrise plot location
- isMale = NA or is_male.txt file (see --counts option)

Set purity & ploidy to `NA` if running unguided.

### `--genome_fa`

Path to genome fasta file.  Expects `.fai` and `.dict` files to be present by appending of extension.

Files extracted from `core_ref*.tar.gz` described [here][cgpwgsrefs].

### `--snp_gc_corr`

Path to file containing SNP loci and GC corrections for ASCAT.

Further information can be found [here][snpgc].

Files extracted from `CNV_SV*.tar.gz` described [here][cgpwgsrefs].

### `--gender_snps`

Path to file representing chrY specific SNP loci.

Absence of these SNPs results in sample being classified as female/XX.  See [here][gensnps] for more details.

Files extracted from `qcGenotype*.tar.gz` described [here][cgpwgsrefs].

### `--counts`

Presence of flag indicates that a previous run has been executed and this is a refinement using pre-generated outputs.

`--pairs` file requires the following under relevant fields:

- reads = `*.count.gz`
- readIdx = `*.count.gz.tbi`
- isMale = `*.is_male.txt`

Use of `purity` and `ploidy` fields is expected.

Only `ascat` step is executed.

## Optional args

### Max job request options

Set the top limit for requested resources for any single job.

| Parameter                             | Description                                                                                                                                                                                                                                                                 | Type     | Default | Required | Hidden |
| ------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ------- | -------- | ------ |
| `max_cpus`                            | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer                                                                            |          |         |          |        |
| e.g. `--max_cpus 1`</small></details> | `integer`                                                                                                                                                                                                                                                                   | 32       |         | True     |        |
| `max_memory`                          | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details> | `string` | 192 GB  |          | True   |
| `max_time`                            | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>        | `string` | 8h      |          | True   |

<!-- refs -->

[cgpwgsrefs]: https://github.com/HealthInnovationEast/cgpwgs-nf?tab=readme-ov-file#reference-files
[gensnps]: https://github.com/cancerit/ascatNgs/wiki#generating-species-specific-gender-loci-file
[snpgc]: https://github.com/cancerit/ascatNgs/wiki#reference-data
