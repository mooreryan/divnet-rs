use crate::config::FitAitchisonConfig;
use crate::conversions;
use crate::conversions::to_composition_matrix;
use crate::matrix;
use crate::matrix::{ewise, Matrix};
use crate::multinomial;
use crate::mvrnorm::mvrnorm;
use blas::dsyrk;
use rand::distributions::Uniform;
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::process;
use std::time::SystemTime;

/// Stores results of fit_aitchison EM iters.
#[derive(Debug)]
pub struct FitAitTmp {
    pub b0: Matrix,
    pub b: Vec<Matrix>,
    pub sigma: Vec<Matrix>,
}

#[derive(Debug)]
pub struct FitAitchsonResult {
    pub beta0: Vec<f64>,
    pub fitted_y: Matrix,
    pub beta: Matrix,
    pub sigma_em: Matrix,
    pub fitted_z: Matrix,
    pub base: usize,
}

// TODO it seems like a pretty high proportion of acceptance...should it be?

/// output should be same dim as fitted_y
fn make_w<R: Rng>(
    rng: &mut R,
    fitted_y: &Matrix,
    sigma: &Matrix,
    sample_sums: &[f64],
    base_taxa: usize,
) -> Matrix {
    let makew_start = SystemTime::now();

    // sigma.ncols() -> NT-1

    let num_samples = fitted_y.ncols();
    let num_taxa = fitted_y.nrows() + 1; // NT-1 w.r.t. W (count matrix)

    assert_eq!(sample_sums.len(), num_samples);

    if sigma.ncols() != fitted_y.nrows() {
        panic!("bad dimensions for sigma and fitted_y");
    }

    let mut new_y = Matrix::zeros(num_taxa - 1, num_samples);

    for sample_idx in 0..num_samples {
        let sample = fitted_y.col(sample_idx);
        // will be 1 x num_taxa

        let norm_vals = mvrnorm(rng, 1, &sample, &sigma);

        assert_eq!((1, num_taxa - 1), norm_vals.dim());

        for (taxa_idx, &x) in norm_vals.row(0).iter().enumerate() {
            new_y.set(taxa_idx, sample_idx, x);
        }
    }

    let compositions = to_composition_matrix(&new_y, base_taxa);
    assert_eq!(compositions.dim(), (num_taxa, num_samples));

    let mut bootstrap_counts = Matrix::zeros(num_taxa, num_samples);

    unsafe {
        for sample_idx in 0..compositions.ncols() {
            //i is also the sample index
            let sample_sum = sample_sums[sample_idx].round() as usize;

            let probs = compositions.col(sample_idx);

            // Each draw specifies a taxa to sample.
            let draws = multinomial::multinomial(rng, 1, sample_sum, &probs);

            for taxa_idx in 0..draws.len() {
                let old_val = bootstrap_counts.get_unchecked(taxa_idx, sample_idx);
                let new_val = old_val + draws[taxa_idx];
                bootstrap_counts.set_unchecked(taxa_idx, sample_idx, new_val);
            }
        }
    }

    log::debug!("`make_w` took {:?}", makew_start.elapsed().unwrap());

    bootstrap_counts
}

pub fn parametric_bootstrap<R: Rng>(
    rng: &mut R,
    fitted_y: &Matrix,
    sigma: &Matrix,
    W: &Matrix,
    X: &Matrix,
    config: &FitAitchisonConfig,
) -> FitAitchsonResult {
    let bootstrap_start = SystemTime::now();

    let mut sample_sums = Vec::new();

    for sample_idx in 0..W.ncols() {
        let sample = W.col(sample_idx);
        sample_sums.push(sample.iter().sum::<f64>());
    }

    //let mw = make_w(&fitted_y, &sigma, sample_sums);

    log::debug!("Starting `make_w`");
    let mw = make_w(rng, &fitted_y, &sigma, &sample_sums, config.base_taxa);

    let res = fit_aitchison(rng, &mw, &X, &config);

    log::debug!(
        "parametric_bootstrap took {:?}",
        bootstrap_start.elapsed().unwrap()
    );

    res
}

fn _diagonal_network(sigma: &Matrix) -> Matrix {
    let inverse: Vec<f64> = sigma.diag().iter().map(|&x| 1. / x).collect();

    Matrix::from_diag(&inverse)
}

fn diagonal_network_vec(sigma: &Matrix) -> Vec<f64> {
    sigma.diag().iter().map(|&x| 1. / x).collect()
}

fn remove_at(v: &[f64], idx: usize) -> Vec<f64> {
    let mut new = Vec::new();
    for i in 0..v.len() {
        if i != idx {
            new.push(v[i])
        }
    }

    new
}

fn _get_Yi<'a>(sample_j: usize, mci: usize, Y: &'a Matrix, Yi_MH: &'a Matrix) -> &'a [f64] {
    if mci == 0 {
        // First mc iter we take the actual Y values for this row.
        Y.col(sample_j)
    } else {
        // The rest of the mc iters, we take the result of the last iteration.
        let last_iteration = mci - 1;

        Yi_MH.col(last_iteration)
    }
}

fn make_Yi_star<R: Rng>(Yi: &[f64], normal: &Normal<f64>, mut rng: &mut R) -> Vec<f64> {
    Yi.iter().map(|&x| x + normal.sample(&mut rng)).collect()
}

// TODO could speed this up by merging the sum and exp.
fn _do_Eq5pt1_old(Yi: &[f64], Yi_star: &[f64], Wi: &[f64]) -> f64 {
    ewise::sum(&Wi)
        * (f64::ln(ewise::sum(&ewise::exp(&Yi)) + 1.)
            - f64::ln(ewise::sum(&ewise::exp(&Yi_star)) + 1.))
}

pub fn do_Eq5pt1(Yi: &[f64], Yi_star: &[f64], Wi: &[f64]) -> f64 {
    let Yi_exp_sum = Yi.iter().map(|&x| f64::exp(x)).sum::<f64>();
    let Yi_star_exp_sum = Yi_star.iter().map(|&x| f64::exp(x)).sum::<f64>();

    Wi.iter().sum::<f64>() * (f64::ln(Yi_exp_sum + 1.) - f64::ln(Yi_star_exp_sum + 1.))
}

// TODO pretty sure this should take Wi_no_base?
pub fn do_Eq5pt2(Yi: &[f64], Yi_star: &[f64], Wi_no_base_taxa: &[f64]) -> f64 {
    let len = Yi.len();
    assert_eq!(len, Yi_star.len());
    assert_eq!(len, Wi_no_base_taxa.len());

    let (w, ystar, y) = (&Wi_no_base_taxa[..len], &Yi_star[..len], &Yi[..len]);

    let mut sum = 0.;

    for i in 0..len {
        sum += w[i] * (ystar[i] - y[i]);
    }

    sum
}

//
//
// This way of doing Eq5pt3 is just like divnet

/*
// R crossprod(X, Y) is t(X) %*% Y
// R: Eq5pt3 <- -0.5 * crossprod((Yi.star - eYi), sigInv) %*% (Yi.star - eYi).
//                                                        ^^^ can just do ewise multiply then sum the result
//
// Make a col vector so we can multiply it.
let Yi_star_minus_eYi =
    Matrix::from_data(1, Yi_star.len(), ewise::sub(&Yi_star, &eYi)).unwrap();

// Yi_star_minus_eYi: length is notus - 1.  sigInv is notus-1 X notus-1
// note that fortran, and Array2 are ColMajor (and so to is BLAS)

let Eq5pt3 = unsafe { Yi_star_minus_eYi.mmuld(&sigInv).unwrap() };
assert_eq!(Eq5pt3.dim(), Yi_star_minus_eYi.dim());
assert_eq!(Eq5pt3.data().len(), Yi_star_minus_eYi.data().len());
let Eq5pt3 = -0.5 * Eq5pt3.ewise_mul(&Yi_star_minus_eYi).sum();

*/
// WARNING...this hack only works because we are are using diagonal networks...in this
// case sigInv is a matrix with elems on the diagonal and zeros everywhere else.  This
// means that Yi_star_minus_eYi and sigInv as vectors and just to elementwise
// multiplication!
fn do_Eq5pt3(Yi_star: &[f64], eYi: &[f64], sigInv_vec: &[f64]) -> f64 {
    assert_eq!(Yi_star.len(), eYi.len());
    let Yi_star_minus_eYi = ewise::sub(&Yi_star, &eYi);
    assert_eq!(Yi_star_minus_eYi.len(), sigInv_vec.len());
    // TODO this (and the orig way can be -inf on the test set.

    -0.5 * (0..Yi_star_minus_eYi.len())
        .map(|i| unsafe {
            let a = Yi_star_minus_eYi.get_unchecked(i);
            a * sigInv_vec.get_unchecked(i) * a
        })
        .sum::<f64>()
}

/*
// R: Eq5pt4 <- -0.5 * crossprod((Yi - eYi), sigInv) %*% (Yi - eYi)
//
// Make a col vector so we can multiply it.
let Yi_minus_eYi = Matrix::from_data(1, Yi.len(), ewise::sub(&Yi, &eYi)).unwrap();

// Yi_star_minus_eYi: length is notus - 1.  sigInv is notus-1 X notus-1
// note that fortran, and Array2 are ColMajor (and so to is BLAS)
let Eq5pt4 = unsafe { Yi_minus_eYi.mmule(&sigInv).unwrap() };

assert_eq!(Eq5pt4.data().len(), Yi_minus_eYi.data().len());

let Eq5pt4 = -0.5 * Eq5pt4.ewise_mul(&Yi_minus_eYi).sum();
*/
fn do_Eq5pt4(Yi: &[f64], eYi: &[f64], sigInv_vec: &[f64]) -> f64 {
    assert_eq!(Yi.len(), eYi.len());
    let Yi_minus_eYi = ewise::sub(&Yi, &eYi);
    assert_eq!(Yi_minus_eYi.len(), sigInv_vec.len());

    -0.5 * (0..Yi_minus_eYi.len())
        .map(|i| unsafe {
            let a = Yi_minus_eYi.get_unchecked(i);
            a * sigInv_vec.get_unchecked(i) * a
        })
        .sum::<f64>()
}

fn get_acceptance(Eq5pt1: f64, Eq5pt2: f64, Eq5pt3: f64, Eq5pt4: f64) -> f64 {
    // R: fullRat <- Eq5pt1 + Eq5pt2 + Eq5pt3 - Eq5pt4
    let full_rat = Eq5pt1 + Eq5pt2 + Eq5pt3 - Eq5pt4;

    let exp_full_rat = f64::exp(full_rat);
    if exp_full_rat < 1. {
        exp_full_rat
    } else {
        1.
    }
}

// todo original DivNet code checks acceptance for NaN.
fn is_accepted<R: Rng>(acceptance: f64, uniform: &Uniform<f64>, mut rng: &mut R) -> bool {
    uniform.sample(&mut rng) < acceptance
}

// TODO lots of LL data write misses.
/// lr_estimates is one row of Yi or Yi_star...ie logratios of NT-1 taxa for a single sample.
fn _update_lr_containers(
    sample_idx: usize,
    mci: usize,
    accepted: bool,
    lr_estimates: &[f64],
    Yi_MH: &mut Matrix,
    mc_iter_logratios: &mut Vec<Matrix>,
) {
    assert_eq!(Yi_MH.ncols(), lr_estimates.len() + 1);

    // The first column is acceptance.
    // TODO don't need to track this.
    if accepted {
        Yi_MH.set(mci, 0, 1.)
    } else {
        Yi_MH.set(mci, 0, 0.)
    }

    for j in 1..Yi_MH.ncols() {
        // Yi_star has one fewer element than Yi_MH.ncols() since the first col of YI_MH
        // is acceptance.
        Yi_MH.set(mci, j, lr_estimates[j - 1]);

        // track result for sigma calculations later
        mc_iter_logratios[mci].set(sample_idx, j - 1, lr_estimates[j - 1]);
    }
}

fn make_mc_iter_lr_container(size: usize, nsamples: usize, ntaxa: usize) -> Vec<Matrix> {
    (0..size)
        .map(|_| Matrix::zeros(ntaxa - 1, nsamples))
        .collect()
}

/// `Y` logratios (NT-1 x NS)
/// `W` counts (NT x NS)
/// `eY` expected logratios (NT-1 x NS)
/// `sigma` (NT-1 x NT-1)
/// `config` configuration values for the function!
fn mcmat<R: Rng>(
    mut rng: &mut R,
    Yi_MH_container: &mut Vec<Matrix>,
    mc_iter_logratios: &mut Vec<Matrix>,
    Y: &Matrix,
    W: &Matrix,
    eY: &Matrix,
    sigma: &Matrix,
    config: &FitAitchisonConfig,
) {
    let ntaxa = W.nrows();
    let nsamples = W.ncols();

    assert_eq!(Y.dim(), (ntaxa - 1, nsamples));
    assert_eq!(eY.dim(), Y.dim());
    assert_eq!(sigma.dim(), (ntaxa - 1, ntaxa - 1));

    // mean 0, standard deviation config.stepsize
    let normal = Normal::new(0.0, config.stepsize).unwrap();
    // Technically, this excludes 1.
    let uniform = Uniform::from(0.0..1.0);

    // Note, this little hack works because the diagonal network is just the entries on the
    // diagonal.  So we can treat it like elementwise operations later.
    // dim is NT-1 x NT-1
    // let sigInv = diagonal_network(&sigma);
    let sigInv_vec = diagonal_network_vec(&sigma);
    assert_eq!(sigInv_vec.len(), ntaxa - 1);

    // TODO don't bother tracking acceptance.
    // This will contain NS Matrices of MCI x NT matrices (the first col of each of
    // these will be acceptance)
    // let mut Yi_MH_container: Vec<Matrix> = Vec::new();

    // // This will contain MCI matrices of NT-1 x NS matrices.  Each one is one MC iteration of
    // // logratio estimates. Need these to update sigma.
    // let mut mc_iter_logratios: Vec<Matrix> =
    //     make_mc_iter_lr_container(config.mc_iter, nsamples, ntaxa);

    log::debug!("Processing samples in MC mat");
    for sample_idx in 0..nsamples {
        log::trace!("Working on sample {} of {}", sample_idx + 1, nsamples);

        // let sample_start_time = SystemTime::now();

        // NT-1 x MCI.  Each taxa has MCI iters of Yi estimates.
        // let mut Yi_MH = Matrix::zeros(ntaxa - 1, config.mc_iter);
        let Yi_MH = &mut Yi_MH_container[sample_idx];
        assert_eq!(Yi_MH.ncols(), config.mc_iter - config.mc_burn);
        assert_eq!(Yi_MH.ncols(), config.mc_burn);

        // todo will need to crash with better msg on weird user input for these!
        assert!(config.mc_iter > 1, "MC iters must be > 0");
        assert!(config.mc_burn > 1, "MC burn must be > 0");
        // im requiring this so this 1/2 mem hack works
        assert_eq!(
            config.mc_iter,
            config.mc_burn * 2,
            "MC iter must be 2 * MC burn "
        );

        for true_mci in 0..config.mc_iter {
            let in_burn_stage = true_mci < config.mc_burn;

            let mci = if in_burn_stage {
                // in the burn stage
                true_mci
            } else {
                // in the keep stage
                true_mci - config.mc_burn
            };

            // Depending on the iteration, either take original Y, or previous Yi, or ___.
            let Yi = if mci == 0 && in_burn_stage {
                // First mc iter we take the actual Y values for this row.
                Y.col(sample_idx)
            } else if mci == 0 {
                // We are NOT in burn stage.  First iter of keep stage.
                //
                // Get the previous Yi_MH.  It is at the end of the storage.
                assert_eq!(
                    config.mc_burn,
                    Yi_MH.ncols(),
                    "mc_burn should equal Yi_MH.cols()"
                );
                let last_iteration = Yi_MH.ncols() - 1;
                Yi_MH.col(last_iteration)
            } else {
                // The rest of the mc iters, we take the result of the last iteration.
                let last_iteration = mci - 1;
                Yi_MH.col(last_iteration)
            };

            let Wi = W.col(sample_idx);
            // This one allocates!
            let Wi_no_base_taxa: Vec<f64> = remove_at(&Wi, config.base_taxa);
            let eYi = eY.col(sample_idx);

            // Offset Yi by draws from a normal distribution.
            // This one allocates!
            let Yi_star: Vec<f64> = make_Yi_star(&Yi, &normal, &mut rng);

            let Eq5pt1 = do_Eq5pt1(&Yi, &Yi_star, &Wi);
            let Eq5pt2 = do_Eq5pt2(&Yi, &Yi_star, &Wi_no_base_taxa);
            let Eq5pt3 = do_Eq5pt3(&Yi_star, &eYi, &sigInv_vec);
            let Eq5pt4 = do_Eq5pt4(&Yi, &eYi, &sigInv_vec);

            let acceptance = get_acceptance(Eq5pt1, Eq5pt2, Eq5pt3, Eq5pt4);
            let accepted = is_accepted(acceptance, &uniform, &mut rng);

            assert_eq!(Yi_MH.nrows(), Yi_star.len());
            assert_eq!(mc_iter_logratios.len(), config.mc_iter - config.mc_burn);
            assert_eq!(mc_iter_logratios[0].dim(), (ntaxa - 1, nsamples));
            if accepted {
                for taxa_idx in 0..Yi_MH.nrows() {
                    let val = Yi_star[taxa_idx];

                    Yi_MH.set(taxa_idx, mci, val);

                    mc_iter_logratios[mci].set(taxa_idx, sample_idx, val);
                }
            } else {
                for taxa_idx in 0..Yi_MH.nrows() {
                    //let val = Yi[taxa_idx];
                    //let val = 10101.;

                    // TODO this extra check is yucky!
                    let val = if mci == 0 && in_burn_stage {
                        // get the original counts...this is truly the first iter.
                        Y.get(taxa_idx, sample_idx)
                    } else if mci == 0 {
                        // We are in the keep stage....the lat iter is the END of the Yi_MH since it
                        // only holds enough space for the real counts.
                        let last_itr = Yi_MH.ncols() - 1;

                        Yi_MH.get(taxa_idx, last_itr)
                    } else {
                        let last_itr = mci - 1;

                        Yi_MH.get(taxa_idx, last_itr)
                    };

                    Yi_MH.set(taxa_idx, mci, val);

                    mc_iter_logratios[mci].set(taxa_idx, sample_idx, val);
                }
            }
        }

        // Yi_MH_container.push(Yi_MH);
    }

    // mc_iter_logratios
}

// TODO make this nicer...it's copied from original code now.
unsafe fn _sigma_sum_function(mc_iter_matrix: &Matrix, eY: &Matrix) -> Matrix {
    let tmp = mc_iter_matrix.ewise_sub(&eY);

    tmp.mmul2(&tmp, matrix::TRANSPOSE, matrix::NO_TRANSPOSE)
        .expect("matrix mult failed ...oops!")
}

fn update_result(ldc: usize, n: usize, c: Vec<f64>) -> Matrix {
    // Now need to make matrix from c.
    let mut result = Matrix::from_data(ldc, n, c).unwrap();

    assert_eq!(result.ncols(), result.nrows());

    // now need to fill in the lower triangular of result.
    for j in 0..result.ncols() {
        for i in 0..result.nrows() {
            if i > j {
                unsafe {
                    result.set_unchecked(i, j, result.get_unchecked(j, i));
                }
            }
        }
    }

    result
}

unsafe fn _sigma_sum_function_dsyrk_ata(mc_iter_matrix: &Matrix, eY: &Matrix) -> Matrix {
    let tmp = mc_iter_matrix.ewise_sub(&eY);

    // use upper part
    let uplo = b'U';

    // means we want A' * A
    let trans = b'T';

    // order of the matrix C (ie the output)
    // A' * A so ncols x nrows X nrows x ncols (of A)
    let n = tmp.ncols();

    // since trans is T, then k is nrows of A
    let k = tmp.nrows();

    // scalar to multiply:  alpha * A' * A
    let alpha = 1.;

    let a = tmp.data();

    // leading dimension is k because trans is 'T'
    let lda = k;

    // scaler: alpha * A' * A + beta * C
    let beta = 1.;

    // this is to hold the output
    let ldc = n;
    let mut c = vec![0.; ldc * n];

    dsyrk(
        uplo, trans, n as i32, k as i32, alpha, &a, lda as i32, beta, &mut c, ldc as i32,
    );

    update_result(ldc, n, c)
}

// This one does AA'
unsafe fn sigma_sum_function_dsyrk_aat(mc_iter_matrix: &Matrix, eY: &Matrix) -> Matrix {
    let tmp = mc_iter_matrix.ewise_sub(&eY);

    // use upper part
    let uplo = b'U';

    // means we want A * A'
    let trans = b'N';

    // order of the matrix C (ie the output)
    // A * A' so nrows x ncols X ncols x nrows (of A ie result is NT-1 x NT-1)
    let n = tmp.nrows();

    // since trans is N, then k is ncols of A
    let k = tmp.ncols();

    // scalar to multiply:  alpha * A' * A
    let alpha = 1.;

    let a = tmp.data();

    // leading dimension is n because trans is 'N'
    let lda = n;

    // scaler: alpha * A' * A + beta * C
    let beta = 1.;

    // this is to hold the output
    let ldc = n;
    let mut c = vec![0.; ldc * n];

    dsyrk(
        uplo, trans, n as i32, k as i32, alpha, &a, lda as i32, beta, &mut c, ldc as i32,
    );

    update_result(ldc, n, c)
}

/// Takes mean value of each taxa across the MCIters.
fn make_Y_new(
    nsamples: usize,
    ntaxa: usize,
    Yi_MH_container: &Vec<Matrix>,
    config: &FitAitchisonConfig,
) -> Matrix {
    // the first column of yi_mh will be acceptance ratio...so let's just ignore that for now
    let mut Y_new = Matrix::zeros(ntaxa - 1, nsamples);
    for (sample_idx, yi_mh) in Yi_MH_container.iter().enumerate() {
        assert_eq!(yi_mh.dim(), (ntaxa - 1, config.mc_iter - config.mc_burn));
        assert_eq!(config.mc_iter - config.mc_burn, config.mc_iter / 2);

        for taxa_idx in 0..yi_mh.nrows() {
            let sample_lr_estimates = yi_mh.row(taxa_idx);

            let sum = sample_lr_estimates //[config.mc_burn..config.mc_iter]
                .iter()
                .sum::<f64>();
            let mean = sum / (config.mc_iter - config.mc_burn) as f64;

            // let mut sum = 0.0;
            // for emi in config.mc_burn..config.mc_iter {
            //     sum += col[emi];
            // }
            //
            // let mean = sum / (config.mc_iter - config.mc_burn) as f64;

            // Y_new has one fewer col since we ignore acceptance!
            Y_new.set(taxa_idx, sample_idx, mean);
        }
    }

    Y_new
}

/// `logratios` should be NT-1 x NS.
///
/// `b0` is the mean taxa LR value across all samples.
fn get_expected_logratios(logratios: &Matrix, b0: &[f64]) -> Matrix {
    // No. taxa present in logratios is one fewer than in the counts!
    let ntaxa = logratios.nrows();
    let nsamples = logratios.ncols();

    assert_eq!(ntaxa, b0.len());

    let mut expected_logratios = Matrix::zeros(ntaxa, nsamples);

    for sample_j in 0..expected_logratios.ncols() {
        for taxa_i in 0..expected_logratios.nrows() {
            let mean_taxa_lr = b0[taxa_i];

            expected_logratios.set(taxa_i, sample_j, mean_taxa_lr);
        }
    }

    expected_logratios
}

// b0 is the mean of each taxa across all samples in the logratio matrix.  It will be indexed with a
// 'j'.
/// `logratios` is NT-1 x NS
pub fn get_b0(logratios: &Matrix) -> Vec<f64> {
    let ntaxa = logratios.nrows() + 1;
    let nsamples = logratios.ncols();

    let mut b0: Vec<f64> = vec![0.; ntaxa - 1];

    for (taxa_i, val) in b0.iter_mut().enumerate() {
        let taxa_mean = logratios.row(taxa_i).iter().sum::<f64>() / nsamples as f64;

        *val = taxa_mean;
    }

    b0
}

/// Like Xc * b = Yc from the paper.
///
/// `centered_covariates` will be NS x NV
/// `centered_logratios` will be NT-1 x NS (so pretty sure you should transpose it!)
///
/// output is b (NV x NT-1)
///
/// TODO switch to doing the transpose in BLAS
pub unsafe fn get_b(centered_covariates: &Matrix, centered_logratios: &Matrix) -> Matrix {
    let b = {
        centered_covariates
            .solve(&centered_logratios.transpose())
            .expect("Failed to solve X_centered * B = Y_centered for B")
    };

    let nvars = centered_covariates.ncols();
    let ntaxa = centered_logratios.nrows() + 1;
    assert_eq!(b.dim(), (nvars, ntaxa - 1));

    b
}

/// Updates b0 with the mean value of each taxa across the samples.
///
/// `Y_new` has NT-1 x NS
/// `b0` has len NT-1;
fn update_b0(b0: &mut [f64], Y_new: &Matrix) {
    assert_eq!(b0.len(), Y_new.nrows());
    for (taxa_idx, elem) in b0.iter_mut().enumerate() {
        let taxa = Y_new.row(taxa_idx);
        let mean = taxa.iter().sum::<f64>() / taxa.len() as f64;

        *elem = mean;
    }
}

fn update_sigma(
    nsamples: usize,
    ntaxa: usize,
    expected_logratios: &Matrix,
    mc_iter_logratios: &Vec<Matrix>,
    sigma: &mut Matrix,
    config: &FitAitchisonConfig,
) {
    assert_eq!(config.mc_iter - config.mc_burn, config.mc_iter / 2);

    // update sigma....here is where we need the NS x NT-1 matrices that are produced by MC
    // iters.  There will be MCiters of them to deal with.
    assert_eq!(mc_iter_logratios.len(), config.mc_iter - config.mc_burn);
    // ignore burnt iters

    assert_eq!(sigma.dim(), (ntaxa - 1, ntaxa - 1));

    log::trace!("`update_sigma` step 1");

    let mut updated_sigma = Matrix::zeros(sigma.nrows(), sigma.ncols());
    for mci in 0..(config.mc_iter - config.mc_burn) {
        // config.mc_burn..mc_iter_logratios.len() {
        let mc_lrs = &mc_iter_logratios[mci];
        assert_eq!(mc_lrs.dim(), (ntaxa - 1, nsamples));
        assert_eq!(expected_logratios.dim(), (ntaxa - 1, nsamples));

        let new_sig = unsafe { sigma_sum_function_dsyrk_aat(&mc_lrs, &expected_logratios) };
        assert_eq!(new_sig.dim(), updated_sigma.dim());

        // add this to some weird temp sigma thingy
        for j in 0..updated_sigma.ncols() {
            for i in 0..updated_sigma.nrows() {
                unsafe {
                    let old_val = updated_sigma.get_unchecked(i, j);

                    updated_sigma.set_unchecked(i, j, old_val + new_sig.get_unchecked(i, j));
                }
            }
        }
    }

    log::trace!("`update_sigma` step 2");

    let mut updated_sigma_mean = Matrix::zeros(sigma.nrows(), sigma.ncols());
    for j in 0..updated_sigma_mean.ncols() {
        for i in 0..updated_sigma_mean.nrows() {
            unsafe {
                let mean = updated_sigma.get_unchecked(i, j)
                    / (nsamples * (config.mc_iter - config.mc_burn)) as f64;
                updated_sigma_mean.set_unchecked(i, j, mean);
            }
        }
    }

    // Actually update sigma.
    *sigma = updated_sigma_mean;
}

pub fn get_sigma_em(num_taxa: usize, fa_tmp: &FitAitTmp, config: &FitAitchisonConfig) -> Matrix {
    // NT-1 x NT-1
    let mut sigma_em = Matrix::zeros(num_taxa - 1, num_taxa - 1);

    for j in 0..sigma_em.ncols() {
        for i in 0..sigma_em.nrows() {
            for emi in (config.em_burn)..(config.em_iter + 1) {
                let mut vals = Vec::new();
                vals.push(fa_tmp.sigma[emi].get(i, j));

                let mean = vals.iter().sum::<f64>() / vals.len() as f64;
                sigma_em.set(i, j, mean);
            }
        }
    }

    sigma_em
}

/// All the OLS stuff is centered so that the mean of LRs across sample is zero.
fn get_fitted_y(
    nsamples: usize,
    ntaxa: usize,
    centered_covariates: &Matrix,
    beta: &Matrix, // OLS result
    beta0: &[f64],
) -> Matrix {
    let nvars = beta.nrows();

    assert_eq!(centered_covariates.nrows(), nsamples);
    assert_eq!(centered_covariates.ncols(), nvars);
    assert_eq!(beta.ncols(), ntaxa - 1);

    // TODO do the transpose in the BLAS code!
    let mut fitted_y = unsafe { centered_covariates.mmulc(&beta).unwrap().transpose() };

    assert_eq!(fitted_y.dim(), (ntaxa - 1, nsamples));
    assert_eq!(beta0.len(), ntaxa - 1);

    for sample_idx in 0..fitted_y.ncols() {
        for taxa_idx in 0..fitted_y.nrows() {
            unsafe {
                let old_val = fitted_y.get_unchecked(taxa_idx, sample_idx);
                let new_val = old_val + beta0[taxa_idx];

                fitted_y.set_unchecked(taxa_idx, sample_idx, new_val);
            }
        }
    }

    fitted_y
}

// TODO this looks wrong....
fn get_beta(
    num_variables: usize,
    num_taxa: usize,
    config: &FitAitchisonConfig,
    fa_tmp: &FitAitTmp,
) -> Matrix {
    let mut beta = Matrix::zeros(num_variables, num_taxa - 1);
    for j in 0..beta.ncols() {
        for i in 0..beta.nrows() {
            for emi in (config.em_burn)..(config.em_iter + 1) {
                let mut vals = Vec::new();
                vals.push(fa_tmp.b[emi].get(i, j));

                let mean = vals.iter().sum::<f64>() / vals.len() as f64;

                unsafe {
                    beta.set_unchecked(i, j, mean);
                }
            }
        }
    }

    beta
}

/// Calculate mean value for all b0 passed the EM burn in.
fn get_beta0(fa_tmp: &FitAitTmp, config: &FitAitchisonConfig) -> Vec<f64> {
    // This is the mean value for all the b0 passed the burn in
    let mut beta0 = Vec::new(); // NT-1

    // this should be emi
    for j in 0..fa_tmp.b0.ncols() {
        // and this should be all taxa estimates for this emi
        let col = fa_tmp.b0.col(j);

        let mut sum = 0.;
        let mut count = 0;

        for (i, &x) in col.iter().enumerate() {
            if i > config.em_burn {
                sum += x;
                count += 1;
            }
        }

        let mean = sum / count as f64;

        beta0.push(mean);
    }

    beta0
}

/// TODO does not validate the values in `config`
pub fn fit_aitchison<R: Rng>(
    rng: &mut R,
    counts: &Matrix,
    covariates: &Matrix,
    config: &FitAitchisonConfig,
) -> FitAitchsonResult {
    log::debug!("Starting `fit_aitchison`");

    let nsamples = counts.ncols();
    let ntaxa = counts.nrows();
    let nvars = covariates.ncols();

    assert_eq!(covariates.dim(), (nsamples, nvars));

    let mut fa_tmp = FitAitTmp {
        b: Vec::new(),
        // TODO should this transpose as well?
        // b0: Matrix::zeros(config.mc_iter + 1, ntaxa - 1),
        b0: Matrix::zeros(config.em_iter + 1, ntaxa - 1),
        sigma: Vec::new(),
    };

    if nsamples != covariates.nrows() {
        log::error!("`nsamples` in count table doesn't match `nsamples` in covariate table");
        process::exit(1);
    }

    let centered_covariates = covariates.center(true);

    let logratios = conversions::to_log_ratios(&counts, config.base_taxa, config.perturbation);

    // b0 is the mean of each taxa across samples in the logratio matrix.  It will be indexed with a
    // 'j'.
    let mut b0 = get_b0(&logratios);
    assert_eq!(b0.len(), ntaxa - 1);

    // NT-1 x NS
    let mut expected_logratios = get_expected_logratios(&logratios, &b0);
    assert_eq!(expected_logratios.dim(), (ntaxa - 1, nsamples));

    // Needed for the OLS solving.
    let centered_logratios = logratios.center(false);

    // R code:  b <- OLS(X_c, Y_p)
    assert_eq!(centered_covariates.dim(), (nsamples, nvars));
    assert_eq!(centered_logratios.dim(), (ntaxa - 1, nsamples));
    let mut b = unsafe { get_b(&centered_covariates, &centered_logratios) };
    assert_eq!(b.dim(), (nvars, ntaxa - 1));

    // R code:  eY <- eY + X_c %*% b
    assert_eq!(centered_covariates.dim(), (nsamples, nvars));
    assert_eq!(b.dim(), (nvars, ntaxa - 1));
    let mut centered_covariates_times_b =
        unsafe { centered_covariates.mmula(&b).unwrap().transpose() };
    assert_eq!(centered_covariates_times_b.dim(), (ntaxa - 1, nsamples));

    assert_eq!(expected_logratios.dim(), centered_covariates_times_b.dim());
    for (i, elem) in expected_logratios.data.iter_mut().enumerate() {
        *elem += centered_covariates_times_b.data[i];
    }

    assert_eq!(logratios.dim(), (ntaxa - 1, nsamples));
    assert_eq!(expected_logratios.dim(), (ntaxa - 1, nsamples));
    // Covariance is a colwise operation.
    let mut sigma = logratios
        .ewise_sub(&expected_logratios)
        .transpose()
        .covariance();
    assert_eq!(sigma.dim(), (ntaxa - 1, ntaxa - 1));

    // TODO transpose fa_tmp.b0
    assert_eq!(b0.len(), ntaxa - 1);
    // store the first of the results
    for (taxa_idx, &x) in b0.iter().enumerate() {
        fa_tmp.b0.set(0, taxa_idx, x);
    }

    fa_tmp.b.push(b.clone());
    fa_tmp.sigma.push(sigma.clone());

    log::debug!("Starting EM iters");
    for em in 0..config.em_iter {
        let em_start = SystemTime::now();

        log::info!("Starting EM iteration {} of {}", em + 1, config.em_iter);

        log::trace!("Running `mcmat`");
        // todo pass in Yi_MH_container?
        // TODO don't bother tracking acceptance.
        // This will contain NS Matrices of MCI x NT matrices (the first col of each of
        // these will be acceptance)
        // let mut Yi_MH_container000: Vec<Matrix> = Vec::new();
        // NT-1 x MCI.  Each taxa has MCI iters of Yi estimates.
        // let mut Yi_MH = Matrix::zeros(ntaxa - 1, config.mc_iter);

        // (NT-1) x nsamples x mc_iter
        let mut Yi_MH_container: Vec<Matrix> = (0..nsamples)
            // todo technically we burn 1/2 of the mc iters so we could do something tricky to only allocate 1/2 and reuse the same part.  we need to keep them tho i think.
            .map(|_i| Matrix::zeros(ntaxa - 1, config.mc_iter - config.mc_burn))
            .collect();

        // This will contain MCI matrices of NT-1 x NS matrices.  Each one is one MC iteration of
        // logratio estimates. Need these to update sigma.
        let mut mc_iter_logratios: Vec<Matrix> =
            make_mc_iter_lr_container(config.mc_iter - config.mc_burn, nsamples, ntaxa);

        // also updates Yi_MH_container
        mcmat(
            rng,
            &mut Yi_MH_container,
            &mut mc_iter_logratios,
            &logratios,
            &counts,
            &expected_logratios,
            &sigma,
            &config,
        );

        log::trace!("Running `make_Y_new`");
        // Estimated logrations from this round of MCiters
        // done 1/2
        let Y_new = make_Y_new(nsamples, ntaxa, &Yi_MH_container, &config);
        assert_eq!(Y_new.dim(), (ntaxa - 1, nsamples));

        log::trace!("Running `update_b0`");
        assert_eq!(b0.len(), ntaxa - 1);
        update_b0(&mut b0, &Y_new);

        log::trace!("Running `update_sigma`");
        // TODO orig code has it this way, but why isn't the updated expected logratios used here?
        // done 1/2
        update_sigma(
            nsamples,
            ntaxa,
            &expected_logratios,
            &mc_iter_logratios,
            &mut sigma,
            &config,
        );

        log::trace!("Setting expected logratios");
        // Note that b0 has been updated by this point.
        for sample_idx in 0..expected_logratios.ncols() {
            for taxa_idx in 0..expected_logratios.nrows() {
                let mean_val_for_this_taxa = b0[taxa_idx];

                expected_logratios.set(taxa_idx, sample_idx, mean_val_for_this_taxa);
            }
        }

        log::trace!("Running `get_b`");
        // TODO this is repeated about 50 lines up...before the em iters start.

        // TODO This one allocates a new b matrix rather than overwriting the old one
        b = unsafe { get_b(&centered_covariates, &Y_new.center(false)) };
        //b = unsafe { get_b(&centered_covariates, &Y_new) }; // todo why not centered here?  idk but it works
        assert_eq!(b.dim(), (nvars, ntaxa - 1));

        log::trace!("Multiplying centered covaiates and b");
        // TODO This one allocates a new b matrix rather than overwriting the old one
        assert_eq!(centered_covariates.dim(), (nsamples, nvars));
        assert_eq!(b.dim(), (nvars, ntaxa - 1));
        // TODO do the transpose in the BLAS code!
        centered_covariates_times_b = unsafe { centered_covariates.mmulb(&b).unwrap().transpose() };
        assert_eq!(centered_covariates_times_b.dim(), (ntaxa - 1, nsamples));

        log::trace!("Adding...");
        // TODO last time expected logratios used b0...why not this time?
        assert_eq!(expected_logratios.dim(), centered_covariates_times_b.dim());
        // elementwise addition here!
        for (i, elem) in expected_logratios.data.iter_mut().enumerate() {
            *elem += centered_covariates_times_b.data[i];
        }

        log::trace!("Setting...");
        // TODO confusing name! in fa_tmp, b0 is EM by taxa
        assert_eq!(b0.len(), ntaxa - 1);
        for (taxa_idx, &x) in b0.iter().enumerate() {
            fa_tmp.b0.set(em + 1, taxa_idx, x);
            //fa_tmp.b0.set(em + 1, taxa_idx, x);
        }

        log::trace!("Pushing...");
        fa_tmp.b.push(b.clone());
        fa_tmp.sigma.push(sigma.clone());

        log::debug!("em iter {} took {:?}", em + 1, em_start.elapsed().unwrap());
    }
    log::trace!("Just finished EM iters");

    // This is the mean value for all the b0 passed the burn in
    let beta0 = get_beta0(&fa_tmp, &config);
    assert_eq!(beta0.len(), ntaxa - 1);

    // NV x NT-1
    let beta = get_beta(nvars, ntaxa, &config, &mut fa_tmp);
    assert_eq!(beta.dim(), (nvars, ntaxa - 1));

    let fitted_y = get_fitted_y(nsamples, ntaxa, &centered_covariates, &beta, &beta0);
    assert_eq!(fitted_y.dim(), (ntaxa - 1, nsamples));

    // // NT-1 x NT-1
    let sigma_em = get_sigma_em(ntaxa, &mut fa_tmp, &config);
    assert_eq!(sigma_em.dim(), (ntaxa - 1, ntaxa - 1));

    let fitted_z = conversions::to_composition_matrix(&fitted_y, config.base_taxa);
    assert_eq!(fitted_z.dim(), (ntaxa, nsamples));

    let thingy = FitAitchsonResult {
        base: config.base_taxa,
        beta0,
        beta,
        fitted_y,
        sigma_em,
        fitted_z,
    };

    log::trace!("Just finished `fit_aitchison`");

    // TODO rename this!
    thingy
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use std::f64;

    const TOL: f64 = 1e-5;

    // All the tests have .transpose() as orignally the data was NS x NT to match DivNet.  But now
    // it processes data as NT x NS.  But I didn't want to change all the examples!

    #[test]
    fn fit_aitchison_gives_reasonable_result() {
        let count_table = Matrix::from_data(
            3,
            3,
            vec![100., 100., 100., 100., 100., 100., 100., 100., 100.],
        )
        .unwrap()
        .transpose();

        let sample_data = Matrix::from_data(3, 1, vec![0., 0., 0.]).unwrap();

        let config = FitAitchisonConfig {
            em_iter: 6,
            em_burn: 3,
            mc_iter: 500,
            mc_burn: 250,
            stepsize: 0.01,
            perturbation: 0.05,
            replicates: 1,
            base_taxa: 0,
        };

        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let result = fit_aitchison(&mut rng, &count_table, &sample_data, &config);

        let expected = vec![0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33, 0.33];

        approx::assert_abs_diff_eq!(result.fitted_z.data(), &expected[..], epsilon = 0.025);
    }

    #[test]
    fn fit_aitchison_gives_reasonable_result2() {
        let count_table = Matrix::from_data(
            4,
            6,
            vec![
                1000., 1000., 0., 0., 1000., 1000., 0., 0., 1000., 1000., 0., 0., 0., 0., 1000.,
                1000., 0., 0., 1000., 1000., 0., 0., 1000., 1000.,
            ],
        )
        .unwrap()
        .transpose();

        let sample_data = Matrix::from_data(4, 1, vec![1., 1., 0., 0.]).unwrap();

        // TODO when you crank up the em iter and mc iter, this actually seems to get less accurate?
        let config = FitAitchisonConfig {
            em_iter: 7,
            em_burn: 3,
            mc_iter: 200,
            mc_burn: 100,
            stepsize: 0.01,
            perturbation: 0.05,
            replicates: 1,
            base_taxa: 0,
        };

        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let result = fit_aitchison(&mut rng, &count_table, &sample_data, &config);

        let expected = vec![
            0.33, 0.33, 0.33, 0., 0., 0., 0.33, 0.33, 0.33, 0., 0., 0., 0., 0., 0., 0.33, 0.33,
            0.33, 0., 0., 0., 0.33, 0.33, 0.33,
        ];

        approx::assert_abs_diff_eq!(result.fitted_z.data(), &expected[..], epsilon = 0.05);
    }

    #[test]
    fn diagonal_network_makes_diagonal_networks() {
        let sigma = Matrix::from_data(2, 2, vec![10., 1., 1., 10.])
            .unwrap()
            .transpose();

        let actual = _diagonal_network(&sigma);

        let expected = Matrix::from_data(2, 2, vec![0.1, 0., 0., 0.1])
            .unwrap()
            .transpose();

        approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn remove_at_returns_a_vec_without_idx() {
        let v = vec![1., 2., 3., 4.];

        approx::assert_relative_eq!(
            &remove_at(&v, 0)[..],
            &vec![2., 3., 4.,][..],
            max_relative = TOL
        );

        approx::assert_relative_eq!(
            &remove_at(&v, 1)[..],
            &vec![1., 3., 4.,][..],
            max_relative = TOL
        );

        approx::assert_relative_eq!(
            &remove_at(&v, 2)[..],
            &vec![1., 2., 4.,][..],
            max_relative = TOL
        );

        approx::assert_relative_eq!(
            &remove_at(&v, 3)[..],
            &vec![1., 2., 3.,][..],
            max_relative = TOL
        );
    }
}
