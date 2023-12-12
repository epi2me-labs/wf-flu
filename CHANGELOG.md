# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.0.1]
###
- Updated docs

## [v1.0.0]
### Changed
- Updated docs
- Cloud readiness

## [v0.0.11]
### Changed
- Output directory structure change
- Adjusted size and positioning of badges to improve layout of typing table
- Nextflow minimum version 23.04.2.

## [v0.0.10]
### Updated
- Update to schema to point to cloud demo

## [v0.0.9]
### Changed
- Update licence to ONT Public License
- Made `--sample_sheet` an optional parameter
- Improved `Typing` table in report to make typing results clearer

## [v0.0.8]
### Fixed
- Poor typing performance after INSaFLU update. Rollback to version used in v0.0.6.

## [v0.0.7]
### Added
- Added workflow-glue
- Added Kit 14 test data
- Abricate results included in the report
### Changed
- CI updates
- Removed ping scripts
- Update check_sample_sheet.py
- `--basecall_cfg` is now used to determine a suitable Medaka model, alternatively provide the name of a model with `--medaka_consensus_model` to override automatic selection
- New fastqingress implementation
### Updated
- INSaFLU v10

## [v0.0.6]
### Fixed
- Updated sample sheet to expect a file

## [v0.0.5]
### Changed
- Updated description in manifest

## [v0.0.4]
### Changed
- Nextflow 22.08 compatibility
- Documentation
### Removed
- conda support

## [v0.0.3]
### Changed
- Heatmap scales according to number of samples.
- Unclassified handled separatley on report.
- Output coverage file.
### Fixed
- Bug where if no unclassified report fails.

## [v0.0.2]
### Added
- Quality score cut-off parameter.
### Changed
- Better processing of metadata.

## [v0.0.1]
### Added
- First release.
