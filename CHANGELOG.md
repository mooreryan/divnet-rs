# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added 

- Makefile for testing the Lee phylum dataset during development.

### Changed

- Updated CHANGELOG.md.
- Updated docs to describe seeding.
- Users must now specify a seed for the random number generator in the config files.
- Switch to the ChaCha20 rng rather than `thread_rng()`

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


[Unreleased]: https://github.com/mooreryan/divnet-rs/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/mooreryan/divnet-rs/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/mooreryan/divnet-rs/releases/tag/v0.1.0
