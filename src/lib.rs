#![allow(non_snake_case)]

// Also need this here for the tests.
extern crate blas;
extern crate lapack;
extern crate openblas_src;

pub mod config;
pub mod conversions;
pub mod fit_aitchison;
pub mod io;
pub mod matrix;
pub mod multinomial;
pub mod mvrnorm;
pub mod opts;
