use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

use crate::config::Config;
use crate::matrix::Matrix;

pub fn read_config(path: &PathBuf) -> Config {
    let f = File::open(path).expect("couldn't open config file");
    let mut f = BufReader::new(f);
    let mut config = String::new();
    f.read_to_string(&mut config).unwrap();

    toml::from_str(&config).expect("failed to deserialize config")
}

/// Reads a OTU table with taxa as rows, samples as columns.
///
/// Output is Matrix with dim == NT x NS.
///
/// TODO this fn does almost no checking of user input!
pub fn read_otu_table(path: &PathBuf) -> (Vec<String>, Vec<String>, Matrix) {
    let f = File::open(path).expect("couldn't open OTU table");
    let f = BufReader::new(f);

    let mut nsamples = 0_usize;

    let mut table_data: Vec<f64> = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();
    let mut taxa_names: Vec<String> = Vec::new();

    for (line_idx, line) in f.lines().enumerate() {
        let line = line.expect("couldn't read line!");

        let ary: Vec<&str> = line.split(',').collect();

        if line_idx == 0 {
            nsamples = ary.iter().len() - 1;

            // First line is the header.
            for (sample_idx, &sample_name) in ary.iter().enumerate() {
                // First row, first column is header for taxa names.  Ignore it.
                if sample_idx > 0 {
                    sample_names.push(sample_name.to_string());
                }
            }
        } else if nsamples != ary.iter().len() - 1 {
            panic!("wrong number of columns for line {}", line_idx + 1);
        } else {
            for (sample_idx, &dat) in ary.iter().enumerate() {
                // First row, first column is header for taxa names.  Ignore it.
                if sample_idx == 0 {
                    taxa_names.push(dat.to_string());
                } else {
                    table_data.push(dat.parse().unwrap());
                }
            }
        }
    }

    let ntaxa = taxa_names.iter().len();

    // TODO ideally, you just read the data in the right way and avoid the transpose.
    // Note that we save the data as NT x NS!
    (
        sample_names,
        taxa_names,
        Matrix::from_data(nsamples, ntaxa, table_data)
            .unwrap()
            .transpose(),
    )
}

/// Reads R model matrix output.
pub fn read_sample_data(path: &PathBuf) -> (Vec<String>, Vec<String>, Matrix) {
    let f = File::open(path).expect("couldn't open sample data table");
    let f = BufReader::new(f);

    let mut nvariables = 0_usize;

    let mut table_data: Vec<f64> = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();
    let mut variable_names: Vec<String> = Vec::new();

    for (line_idx, line) in f.lines().enumerate() {
        let line = line.expect("couldn't read line!");

        let ary: Vec<&str> = line.split(',').collect();

        if line_idx == 0 {
            nvariables = ary.iter().len() - 1;

            // First line is the header.
            for (sample_idx, &sample_name) in ary.iter().enumerate() {
                // First row, first column is header for taxa names.  Ignore it.
                if sample_idx > 0 {
                    variable_names.push(sample_name.to_string());
                }
            }
        } else if nvariables != ary.iter().len() - 1 {
            panic!("wrong number of columns for line {}", line_idx + 1);
        } else {
            for (sample_idx, &dat) in ary.iter().enumerate() {
                // First row, first column is header for taxa names.  Ignore it.
                if sample_idx == 0 {
                    sample_names.push(dat.to_string());
                } else {
                    table_data.push(dat.parse().unwrap());
                }
            }
        }
    }

    let nsamples = sample_names.iter().len();

    (
        sample_names,
        variable_names,
        Matrix::from_data(nvariables, nsamples, table_data)
            .unwrap()
            .transpose(),
    )
}

/*
#[cfg(test)]
mod tests {
    use super::*;
    use std::f64;

    #[macro_use]
    use approx;

    const TOL: f64 = 1e-5;

    #[test]
    fn read_otu_table_parses_an_otu_table() {
        approx::assert_relative_eq!(1., 1., max_relative = TOL);
    }
}
*/
