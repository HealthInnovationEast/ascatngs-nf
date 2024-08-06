# ascat-nf

Nextflow pipeline for using [ascatNgs]

## Docs

- [Usage](docs/Usage.md)
- [Local Testing](docs/Local_Testing.md)

## Version mapping

| ascatNgs | ascatngs-nf |
| -------: | ----------: |
|    4.5.0 |       1.0.0 |

See [ascatNgs changes][ascatngs-changes] or [ascatNgs releases][ascatngs-releases]

## Development and Release process

Please use HubFlow methodology to prepare features, releases and hotfixes.

### Finalising a release

1. Update version in [`nextflow.config`](nextflow.config)
1. Update [Version mapping](#version-mapping) table
1. Complete [`CHANGES.md`](CHANGES.md)

<!-- refs -->

[ascatngs]: https://github.com/cancerit/ascatNgs
[ascatngs-changes]: https://github.com/cancerit/ascatNgs/blob/master/CHANGES.md
[ascatngs-releases]: https://github.com/cancerit/ascatNgs/releases
