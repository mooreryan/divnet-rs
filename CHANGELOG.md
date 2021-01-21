# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Changed

- Move documentation from the `./docs` directory into `./book/src`.  This sets up using `mdBook` for generating nice online documentation.

## [0.2.0] -- 2020-01-19

### Added

- Makefile for testing the Lee phylum dataset during development.
- Documentation for:
    - Properly setting MC iter and MC burn options
    - Setting OpenBLAS threads
    - Getting sample order correct in the sample data file

### Changed

- Updated CHANGELOG.md.
- Updated docs to describe seeding.
- Users must now specify a seed for the random number generator in the config files.
- Switched to the ChaCha20 rng rather than `thread_rng()`
- Changed the way the MC iterations are stored.  Now `divnet-rs` uses 1/3 less RAM than before!
- Updated the `rand` crate dependencies.

## [0.1.1] -- 2020-01-16

### Added 

- Better installation instructions
- An full worked example using the Lee dataset from DivNet
- More documentation

### Changed

- When compiling, `divnet-rs` now uses the system installed OpenBLAS library rather that using the one bundled with `openblas-src`.

### Fixed

- The eigenvector test failed on certain machines because the sign of the eigenvectors is arbitrary.  This has been fixed.

## [0.1.0] -- 2021-01-14

Initial commit!


[Unreleased]: https://github.com/mooreryan/divnet-rs/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/mooreryan/divnet-rs/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/mooreryan/divnet-rs/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/mooreryan/divnet-rs/releases/tag/v0.1.0
