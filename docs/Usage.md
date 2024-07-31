# Usage <!-- omit in toc -->

!! This is a template document, please update when creating your workflow. !!

All links here are pinned to the version of XXXX that this nextflow has been created for.

- [Required args](#required-args)
  - [`--some_required_param`](#--some_required_param)
- [Optional args](#optional-args)
  - [`--some_optional_param`](#--some_optional_param)
- [Resource args](#resource-args)
  - [`--some_resource_param`](#--some_resource_param)

## Required args

Workflow would fail if these are not defined

### `--inputs`

CSV file with header of:

```
groupId,sampleId,type,protocol,platform,reads,readIdx,purity,ploidy
```

- groupId = Links tumour & normal samples
- sampleId = Name of the sample
- type = case/control
- protocol = e.g. WGS, WES
- platform = e.g. ILLUMINA
- reads = BAM/CRAM file
- readIdx = BAM/CRAM index file
- purity = NA or Purity (rho) setting for manual setting of sunrise plot location
- ploidy = NA or Ploidy (psi) setting for manual setting of sunrise plot location

Set purity & ploidy to `NA` if running unguided.

## Optional args

Expect these to have a default behaviour if not defined

### `--some_optional_param`

Description & conditions

## Resource args

These change compute resources, defaults will be present

### `--some_resource_param`

Description & conditions
