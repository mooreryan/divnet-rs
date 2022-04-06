use crate::matrix::Matrix;
use rand::Rng;
use rand_distr::StandardNormal;
use std::process;

// ~18% C.E.
/// Result will be mu.len() X num_observations Array2.  Remember it's stored in Col major, so it may
/// or may not be optimal.
pub fn mvrnorm<R: Rng>(rng: &mut R, num_observations: usize, mu: &[f64], sigma: &Matrix) -> Matrix {
    let num_variables = mu.len();

    if sigma.dim() != (num_variables, num_variables) {
        log::error!(
            "Incompatible arguments:  sigma.dim() != (mu.len(), mu.len()). \
        sigma.dim(): {:?}, mu.len() (aka num_variables): {}",
            sigma.dim(),
            mu.len()
        );
        process::exit(1);
    }

    let normal_values = Matrix::from_data(
        num_variables,
        num_observations,
        standard_normal(rng, num_variables * num_observations),
    )
    .unwrap();

    let mut result = unsafe {
        match sigma.chol() {
            Ok(chol_lower) => chol_lower
                .mmul(&normal_values)
                .expect("mmul failed at chol_lower * normal_values"),
            Err(_chol_err) => {
                match sigma.eigen() {
                    Ok(eigen_result) => {
                        // TODO do the thing with the multiplying (see mvrnorm)
                        let vals: Vec<f64> = eigen_result
                            .values
                            .iter()
                            .map(|&x| {
                                // Some of the eigenvalues are very likely to be zero since the
                                // Cholesky decomp failed!
                                if x < 0. {
                                    0.
                                } else {
                                    f64::sqrt(x)
                                }
                            })
                            .collect();

                        let vals = Matrix::from_diag(&vals);
                        let tmp = eigen_result.vectors.mmul(&vals).unwrap();

                        tmp.mmul(&normal_values).unwrap()
                    }
                    Err(_eigen_err) => {
                        panic!("the spectral decomp also failed :(");
                    }
                }
            }
        }
    };

    // j is also the variable index
    for j in 0..result.ncols() {
        for i in 0..result.nrows() {
            let mu = mu[i];
            let new_val = result.get(i, j) + mu;

            result.set(i, j, new_val);
        }
    }

    result.transpose()
}

fn standard_normal<R: Rng>(rng: &mut R, n: usize) -> Vec<f64> {
    (0..n).map(|_i| rng.sample(StandardNormal)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx;
    use rand::SeedableRng;
    use rand_pcg::Pcg64;
    use std::f64;

    // TODO this test could give some false positives depending on the random numbers generated!
    #[test]
    fn mvrnorm_gives_correct_mu_and_sigma() {
        let mu = vec![0.3427267, 5.2832232];
        let sigma =
            Matrix::from_data(2, 2, vec![0.5165287, -0.2067872, -0.2067872, 2.5308647]).unwrap();

        let mut rng = Pcg64::seed_from_u64(1);

        // Need a lot of samples to get near the originals.
        let result = mvrnorm(&mut rng, 10000, &mu, &sigma);

        let colmeans: Vec<f64> = (0..result.ncols())
            .map(|j| result.col(j).iter().sum::<f64>() / result.nrows() as f64)
            .collect();

        // Means are approx equal to starting means.
        approx::assert_abs_diff_eq!(&colmeans[..], &mu[..], epsilon = 0.25);

        let covariance = result.covariance();
        assert_eq!(covariance.dim(), sigma.dim());
        approx::assert_abs_diff_eq!(covariance.data(), sigma.data(), epsilon = 0.25);
    }
}
