# template-nf

This is a template project to be used when creating nextflow for CYNAPSE projects.

## Using the template

Nextflow items to be changed:

- `main.nf`
  - Usage block URL updates
  - Search for `Header log info`, add params
  - Define processes and tie together in workflow
- `nextflow.config`
  - URLs
  - Versions
  - Container references (in process resource section)
  - Addition of params & defaults
- `docs/Usage.md`
- `conf/test.config`
  - As appropriate for your example exec
