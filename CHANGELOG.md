# Changelog

## Unreleased

- Switched from `ChaCha20` to `Pcg64`.  This will change the output as compared to earlier versions even with the same seed.
- Changed lots of little implementation details.
  - On the Lee dataset with default settings, as compared version `0.2.1` it is...
  - About 3x faster 
  - And uses 60% of the memory 

## 0.2.1 (2021-01-22)

- Reduced memory usage.
- Fixed some Clippy warnings.

## 0.2.0 (2021-01-19)

- Changed the way the MC iterations are stored.  Now `divnet-rs` uses 1/3 less RAM than before!
- Updated the `rand` crate dependencies.
- Makefile for testing the Lee phylum dataset during development.
- Documentation for:
    - Properly setting MC iter and MC burn options
    - Setting OpenBLAS threads
    - Getting sample order correct in the sample data file
    - Seeding the random number generator

## 0.1.1 (2021-01-16)

- Better installation instructions
- Add a full worked example using the Lee dataset from DivNet
- Add documentation
- When compiling, `divnet-rs` now uses the system installed OpenBLAS library rather that using the one bundled with `openblas-src`.
- The eigenvector test failed on certain machines because the sign of the eigenvectors is arbitrary.  This has been fixed.

## 0.1.0 (2021-01-14)

Initial commit!
