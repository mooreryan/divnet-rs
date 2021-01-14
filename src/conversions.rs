use crate::matrix::Matrix;
use std::cmp::Ordering;

//   out <- matrix(nrow = nn, ncol = qq + 1)
//   out[,-base] <- exp(Y)
//   out[,base] <- 1
//   out / rowSums(out)
/// `logratios` will be NT-1 x NS
pub fn to_composition_matrix(logratios: &Matrix, base: usize) -> Matrix {
    // Composition matrix has one more taxa than logratios.

    let nsamples = logratios.ncols();
    let ntaxa = logratios.nrows() + 1;

    let mut result = Matrix::zeros(ntaxa, nsamples);
    let lr_exp = logratios.ewise_exp();

    // j is samples
    for sample_idx in 0..result.ncols() {
        // i is taxa
        for taxa_idx in 0..result.nrows() {
            match taxa_idx.cmp(&base) {
                // The base taxa gets set to 1 --TODO: double check that this is right!
                Ordering::Equal => result.set(taxa_idx, sample_idx, 1.),
                Ordering::Greater => {
                    // We're currently passed the base taxa.  Need to back the taxa index by one
                    // since
                    // we skip the base taxa.
                    result.set(taxa_idx, sample_idx, lr_exp.get(taxa_idx - 1, sample_idx))
                }
                Ordering::Less => {
                    result.set(taxa_idx, sample_idx, lr_exp.get(taxa_idx, sample_idx))
                }
            }
        }
    }

    // Samples are columns.
    let mut sample_sums = Vec::new();
    for j in 0..result.ncols() {
        let sample = result.col(j);

        sample_sums.push(sample.iter().sum::<f64>());
    }

    assert_eq!(sample_sums.len(), nsamples);

    for j in 0..result.ncols() {
        for i in 0..result.nrows() {
            result.set(i, j, result.get(i, j) / sample_sums[j]);
        }
    }

    result
}

/// `counts` is a NT x NS matrix.
pub fn to_log_ratios(counts: &Matrix, base: usize, perturbation: f64) -> Matrix {
    let ntaxa = counts.nrows();
    let nsamples = counts.ncols();

    // Each element in this is the count of the base taxon in sample j.
    let mut base_taxa: Vec<f64> = counts.row(base);

    assert_eq!(base_taxa.len(), nsamples);

    // Replace zero counts with the perturbation.
    for count in base_taxa.iter_mut() {
        if *count == 0. {
            *count = perturbation;
        }
    }

    let mut logratios = Matrix::zeros(ntaxa - 1, nsamples);

    for sample_j in 0..counts.ncols() {
        for taxa_i in 0..counts.nrows() {
            if taxa_i != base {
                let count_ij = counts.get(taxa_i, sample_j);

                let logratio = if count_ij == 0. {
                    // Replace 0 value with perturbation.
                    f64::ln(perturbation / base_taxa[sample_j])
                } else {
                    f64::ln(count_ij / base_taxa[sample_j])
                };

                if taxa_i < base {
                    logratios.set(taxa_i, sample_j, logratio);
                } else {
                    // We are passed the base now.  Move taxa index back one.
                    logratios.set(taxa_i - 1, sample_j, logratio);
                }
            }
        }
    }

    logratios
}

// TODO: add test case when base taxa is not 0 (index).
#[cfg(test)]
mod tests {
    use super::*;
    use approx;
    use std::f64;

    const TOL: f64 = 1e-5;
    const PERTERBATION: f64 = 0.05;
    const BASE_TAXA: usize = 0;

    #[test]
    fn to_log_ratios_works() {
        let counts = Matrix::from_data(
            3,
            4,
            vec![
                1., 2., 3., 10., 20., 30., 100., 200., 300., 1000., 2000., 3000.,
            ],
        )
        .unwrap()
        .transpose();

        let logratios = to_log_ratios(&counts, BASE_TAXA, PERTERBATION);

        // This values are from R....
        //
        // W <- matrix(c(1,2,3, 10,20,30, 100,200,300, 1000,2000,3000), nrow=3)
        //
        // base <- 1
        //
        // Y <- log(W[, -base] / W[, base])
        let expected = Matrix::from_data(
            3,
            3,
            vec![
                2.302585, 2.302585, 2.302585, 4.60517, 4.60517, 4.60517, 6.907755, 6.907755,
                6.907755,
            ],
        )
        .unwrap()
        .transpose();

        approx::assert_abs_diff_eq!(logratios.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn to_log_ratios_handles_zero_counts_in_base_taxa() {
        let counts = Matrix::from_data(
            3,
            4,
            vec![
                0., 2., 3., 10., 20., 30., 100., 200., 300., 1000., 2000., 3000.,
            ],
        )
        .unwrap()
        .transpose();

        let logratios = to_log_ratios(&counts, BASE_TAXA, 1.);

        // This values are from R....
        //
        // W <- matrix(c(1,2,3, 10,20,30, 100,200,300, 1000,2000,3000), nrow=3)
        //
        // base <- 1
        //
        // Y <- log(W[, -base] / W[, base])
        let expected = Matrix::from_data(
            3,
            3,
            vec![
                2.302585, 2.302585, 2.302585, 4.60517, 4.60517, 4.60517, 6.907755, 6.907755,
                6.907755,
            ],
        )
        .unwrap()
        .transpose();

        approx::assert_abs_diff_eq!(logratios.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn to_log_ratios_handles_zero_counts_in_non_base_taxa() {
        let counts = Matrix::from_data(
            3,
            4,
            vec![
                1., 2., 3., 10., 20., 30., 100., 200., 300., 1000., 0., 3000.,
            ],
        )
        .unwrap()
        .transpose();

        let logratios = to_log_ratios(&counts, BASE_TAXA, 2000.);

        // This values are from R....
        //
        // W <- matrix(c(1,2,3, 10,20,30, 100,200,300, 1000,2000,3000), nrow=3)
        //
        // base <- 1
        //
        // Y <- log(W[, -base] / W[, base])
        let expected = Matrix::from_data(
            3,
            3,
            vec![
                2.302585, 2.302585, 2.302585, 4.60517, 4.60517, 4.60517, 6.907755, 6.907755,
                6.907755,
            ],
        )
        .unwrap()
        .transpose();

        approx::assert_abs_diff_eq!(logratios.data(), expected.data(), epsilon = TOL);
    }

    mod output_matches_r {
        use super::*;

        // The following matrices come from R and DivNet!
        fn make_count_table() -> Matrix {
            Matrix::from_data(2, 4, vec![10., 9., 15., 0., 3., 2., 0., 7.])
                .unwrap()
                .transpose()
        }

        fn _make_freq_table() -> Matrix {
            Matrix::from_data(
                2,
                4,
                vec![
                    0.3571429, 0.5, 0.5357143, 0., 0.1071429, 0.1111111, 0., 0.3888889,
                ],
            )
            .unwrap()
            .transpose()
        }

        fn make_logratio_table() -> Matrix {
            Matrix::from_data(
                2,
                3,
                vec![
                    0.4054651, -5.1929569, -1.203973, -1.504077, -5.2983174, -0.2513144,
                ],
            )
            .unwrap()
            .transpose()
        }

        fn make_composition_table() -> Matrix {
            Matrix::from_data(
                2,
                4,
                vec![
                    0.3565062,
                    0.4986150,
                    0.534759358,
                    0.002770083,
                    0.1069519,
                    0.1108033,
                    0.001782531,
                    0.387811634,
                ],
            )
            .unwrap()
            .transpose()
        }

        #[test]
        fn to_log_ratios_matches_r() {
            let count_table = make_count_table();
            let logratios = make_logratio_table();

            let actual = to_log_ratios(&count_table, BASE_TAXA, PERTERBATION);

            assert_eq!(logratios.dim(), actual.dim());
            approx::assert_abs_diff_eq!(logratios.data(), actual.data(), epsilon = TOL);
        }

        #[test]
        fn to_composition_matrix_matches_r() {
            let logratios = make_logratio_table();
            let compositions = make_composition_table();

            let actual = to_composition_matrix(&logratios, BASE_TAXA);

            assert_eq!(compositions.dim(), actual.dim());
            approx::assert_abs_diff_eq!(compositions.data(), actual.data(), epsilon = TOL);
        }

        #[test]
        fn you_get_the_correct_compositions_from_counts() {
            let counts = make_count_table();

            let logratios = to_log_ratios(&counts, BASE_TAXA, PERTERBATION);
            let compositions = to_composition_matrix(&logratios, BASE_TAXA);

            let expected = make_composition_table();

            assert_eq!(compositions.dim(), expected.dim());
            approx::assert_abs_diff_eq!(compositions.data(), compositions.data(), epsilon = TOL);
        }
    }
}
