# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.11.7] - TBD
- Genomes are now saved in batches to the workspace
- Unusable `export_genome_features_protein_to_fasta` function was removed
- The `genbanks_to_genomes` method was added to allow users to upload multiple
genome objects at once

## [0.11.6] - 2022-01-29
Fixed performance issue with FastaGFF impacting metagenome uploads

## [0.11.5] - 2021-3-16
Fixed issue with empty GFF download file PTV-1637

## [0.11.4] - 2021-2-2
Fixed workspace and handle permission issue

## [0.11.3] - 2020-11-24

### Fixed/Changed
  Rename CHANGELOG.md to RELEASE_NOTES.md 


## [0.11.2] - 2020-10-30

### Fixed/Changed
  Fixed false ontology_ref requirment


## [0.11.1] - 2020-05-27

### Fixed/Changed
- behavior of the Genbank upload has been fixed to accurately load Ensembl data from users.
- Duplicate ID's are now accpeted in genbank files with the string '.n' appended to duplicates. n = "the number of times this id has been seen previously"


## [0.11.0] - 2020-02-18

### Added

- The `update_taxon_assignments` method was added to allow users to create,
  replace, or remove taxon assignments for a Genome object.

## [0.10.4] - 2020-02-20

### Changed
- exporter function for Metagenome fasta + gff added.

## [0.10.3] - 2020-02-20

### Fixed

- ws_obj_gff_to_genome and ws_obj_gff_to_metagenome bug on certain input types resolved

## [0.10.2] - 2020-02-11

### Fixed

- keyError: "warnings" fix for Metagenome upload

## [0.10.1] - 2020-02-10

### Fixed

- ws_obj_gff_to_metagenome was erroneously outputting non-metagenome filetypes when the "is_metagenome" parameter was not specified.
- protein sequences were still sometimes saved for some metagenome objects within the parent gene of a cds. protein_translation no longer saved in metagenome.

## [0.10.0] - 2020-02-05

### Changed

- Adding two new functions: ws_obj_gff_to_metagenome, ws_obj_gff_to_genome

## [0.9.0] - 2019-10-24

### Changed

- Use the Relation Engine API for finding and populating taxonomy data
- Use the NCBI taxonomy ID as the primary way to assign taxonomy data instead of sciname

## [0.8.10] - 2018-11

### Changed

- Use WSLargeDataIO for saving and pulling genomes
- Expand and refactor import of ontology ids

### Fixed

- Fix a few genome warnings

## [0.8.9] - 2018-06-04

### Fixed

- Correct genbank download of old genomes with contig sets and long contig IDs
- Make reference path though genome for assembly access to prevent permissions issues

## [0.8.8] - 2018-05-30

### Changed

- Use searchapi2 in production & appdev

## [0.8.7] - 2018-05-28

### Fixed

- Fix go term splitting for GFF uploads

### Added

- Add optional "upgrade" parameter to save_genome

## [0.8.6] - 2018-05-15

### Changed

- Genome object refactor. Baseline for this log
