use blas::dgemm;
use lapack::{dgels, dpotrf, dsyevr};
use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Sub};

pub mod ewise;

use std::fmt;

// #[macro_use]
use approx;

pub const TRANSPOSE: u8 = b'T';
pub const NO_TRANSPOSE: u8 = b'N';

pub const TOL: f64 = 1e-6;

#[derive(Debug, Clone)]
pub struct EigenOutput {
    pub values: Vec<f64>,
    pub vectors: Matrix,
}

#[derive(Debug, Clone)]
/// A matrix with ColMajor data storage
pub struct Matrix {
    /// `data` is owned by the struct.  Size isn't known at compile time.
    pub data: Vec<f64>,
    pub nrows: usize,
    pub ncols: usize,
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "A {} x {} Matrix", self.nrows, self.ncols).unwrap();

        for i in 0..self.nrows {
            for j in 0..self.ncols {
                write!(f, "{:10.3} ", self.get(i, j)).unwrap();
            }
            writeln!(f).unwrap();
        }

        writeln!(f)
    }
}

// Initializing matrices
impl Matrix {
    /// Initialized data with all zeros.
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        Matrix {
            data: vec![0.; nrows * ncols],
            nrows,
            ncols,
        }
    }

    /// Uses `data` to init the Array2.  Set it up in col major form!
    pub fn from_data(nrows: usize, ncols: usize, data: Vec<f64>) -> Result<Self, String> {
        if data.len() == nrows * ncols {
            Ok(Matrix { data, nrows, ncols })
        } else {
            Err("data.len() != nrows * ncols".to_string())
        }
    }
}

// Getting info about the matrix
impl Matrix {
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }

    pub fn ncols(&self) -> usize {
        self.ncols
    }

    pub fn size(&self) -> usize {
        self.nrows * self.ncols
    }

    pub fn dim(&self) -> (usize, usize) {
        (self.nrows, self.ncols)
    }
}

// Equality
impl Matrix {
    pub fn eq(&self, other: &Self) -> bool {
        self.dim() == other.dim() && self.data == other.data
    }

    pub fn approx_eq(&self, other: &Self, tol: f64) -> bool {
        self.dim() == other.dim()
            && approx::abs_diff_eq!(&self.data[..], &other.data[..], epsilon = tol)
    }
}

// Functions for summarizing aspects of the matrix
impl Matrix {
    pub fn sum(&self) -> f64 {
        self.data.iter().sum::<f64>()
    }

    pub fn row_sums(&self) -> Vec<f64> {
        (0..self.nrows)
            .map(|i| self.row(i).iter().sum::<f64>())
            .collect()
    }

    pub fn col_sums(&self) -> Vec<f64> {
        (0..self.ncols)
            .map(|i| self.col(i).iter().sum::<f64>())
            .collect()
    }
}

macro_rules! define_unary_ops {
    ($($fn:ident => $mod:ident::$op:ident),*) => {$(
        pub fn $fn(&self) -> Self {
            self.ewise_map($mod::$op)
        }
    )*}
}

// Element wise operations for single operands
impl Matrix {
    fn ewise_map(&self, f: fn(f64) -> f64) -> Self {
        Matrix {
            data: ewise::map(f, &self.data),
            nrows: self.nrows,
            ncols: self.ncols,
        }
    }

    define_unary_ops!(ewise_exp => f64::exp, ewise_ln => f64::ln);
}

macro_rules! define_binary_ops {
    ($($fn:ident => $mod:ident::$op:ident),*) => {$(
        pub fn $fn(&self, other: &Self) -> Self {
            self.ewise_zipmap(other, $mod::$op)
        }
    )*}
}

// Element wise operations for two operands.
impl Matrix {
    /// Will panic if the dimensions of Matrices are not the same.
    fn ewise_zipmap(&self, other: &Self, f: fn(f64, f64) -> f64) -> Self {
        assert_eq!(self.dim(), other.dim());

        Matrix {
            data: ewise::zipmap(f, &self.data, &other.data),
            nrows: self.nrows,
            ncols: self.ncols,
        }
    }

    define_binary_ops!(ewise_add => Add::add, ewise_div => Div::div, ewise_mul => Mul::mul, ewise_sub => Sub::sub);
}

// TODO: consider switching checked array access to unchecked
// Get and set elements, columns, rows, etc.
impl Matrix {
    /// Get element ij of self.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        // This is ColMajor access!
        self.data[i + j * self.nrows]
    }

    /// No bounds checking on this one!
    pub unsafe fn get_unchecked(&self, i: usize, j: usize) -> f64 {
        *self.data.get_unchecked(i + j * self.nrows)
    }

    /// Set the value at `i`, `j` to `new_val`.
    pub fn set(&mut self, i: usize, j: usize, new_val: f64) {
        self.data[i + j * self.nrows] = new_val;
    }

    /// # Safety
    ///
    /// This function uses unchecked array access.  Caller must verify that i and j are inbounds.   
    pub unsafe fn set_unchecked(&mut self, i: usize, j: usize, new_val: f64) {
        *self.data.get_unchecked_mut(i + j * self.nrows) = new_val;
    }

    /// Get a "column" of the array as a slice.
    pub fn col(&self, j: usize) -> &[f64] {
        let index = j * self.nrows;

        &self.data[index..index + self.nrows]
    }

    /// Get a "column" of the array as a slice.
    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        let index = j * self.nrows;

        &mut self.data[index..index + self.nrows]
    }

    pub fn multiple_cols(&self, start_j: usize, num_cols: usize) -> &[f64] {
        // TODO: check num_cols & give good error rather than let it panic on array access
        let index = start_j * self.nrows;

        &self.data[index..index + (self.nrows * num_cols)]
    }

    /// Rows are non-contiguous so it returns a new array.
    ///
    /// TODO: provide a no-copy view of the row instead.
    pub fn row(&self, i: usize) -> Vec<f64> {
        // let mut v: Vec<f64> = Vec::new();

        let mut v = vec![0.; self.ncols];

        assert!(i < self.nrows);

        // this way ~50% faster than .iter_mut().enumerate()
        for j in 0..self.ncols {
            unsafe {
                // let mut old_val = v.get_unchecked_mut(j);

                *v.get_unchecked_mut(j) = self.get_unchecked(i, j);
            }
        }

        v
    }
}

impl Matrix {
    // This weird set of duplicated functions with slightly different names are for profiling.  The
    // different versions are used in different spots to be able to figure out which exact matrix
    // multiplation is actually taking up the most runtime.
    //
    // TODO clean this up...remove in production...need a better way to do this!

    // Wrapper when neither is transposed.
    pub unsafe fn mmul(&self, other: &Self) -> Result<Matrix, String> {
        self.mmul2(other, NO_TRANSPOSE, NO_TRANSPOSE)
    }
    pub unsafe fn mmula(&self, other: &Self) -> Result<Matrix, String> {
        self.mmul2_apple(other, NO_TRANSPOSE, NO_TRANSPOSE)
    }
    pub unsafe fn mmulb(&self, other: &Self) -> Result<Matrix, String> {
        self.mmul2_bonky(other, NO_TRANSPOSE, NO_TRANSPOSE)
    }
    pub unsafe fn mmulc(&self, other: &Self) -> Result<Matrix, String> {
        self.mmul2_carlisle(other, NO_TRANSPOSE, NO_TRANSPOSE)
    }
    pub unsafe fn mmuld(&self, other: &Self) -> Result<Matrix, String> {
        self.mmul2_dawn(other, NO_TRANSPOSE, NO_TRANSPOSE)
    }
    pub unsafe fn mmule(&self, other: &Self) -> Result<Matrix, String> {
        self.mmul2_edgy(other, NO_TRANSPOSE, NO_TRANSPOSE)
    }
    pub unsafe fn mmulf(&self, other: &Self) -> Result<Matrix, String> {
        self.mmul2_funzilla(other, NO_TRANSPOSE, NO_TRANSPOSE)
    }

    // unsafe {
    //     dgemm(b'N', b'N', m, n, k, 1.0, &a, m, &b, k, 1.0, &mut c, m);
    // }
    //
    // a <- matrix(c(1., 10., 2., 20., 3., 30.), nrow = 2, ncol = 3)
    // b <- matrix(c(1., 10., 100., 2., 20., 200., 3., 30., 300., 4., 40., 400.), nrow = 3, ncol = 4)
    //
    // a %*% b
    //
    // TODO is there a return code to check?
    /// This is matrix multiplication.
    pub unsafe fn mmul2(
        &self,
        other: &Matrix,
        transpose_self: u8,
        transpose_other: u8,
    ) -> Result<Matrix, String> {
        // TODO check if the matrices can be multiplied
        let transpose_a = transpose_self;
        let transpose_b = transpose_other;
        let mat_a = &self.data;
        let mat_b = &other.data;
        let alpha = 1.;
        let beta = 1.;

        // m is nrows of Op(A) and nrows matrix C
        let m = if transpose_self == TRANSPOSE {
            // We're transposing self aka the a matrix
            self.ncols as i32
        } else {
            self.nrows as i32
        };

        // n is ncols of Op(b) and ncols of matrix c
        let n = if transpose_other == TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        // k is ncols of op(A), and nrows of op(B)
        let k = if transpose_self == TRANSPOSE && transpose_other == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_self == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_other == TRANSPOSE {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        };

        let c_nrows = m;
        let c_ncols = n;

        let mut mat_c = vec![0.; (c_nrows * c_ncols) as usize];

        let lda = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        let ldb = if transpose_other == NO_TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        let ldc = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        // DGEMM: C = alpha * op(A) * op(B) + beta * C.  Note: C gets overwritten with this.

        dgemm(
            transpose_a,
            transpose_b,
            m,
            n,
            k,
            alpha,
            mat_a,
            lda,
            mat_b,
            ldb,
            beta,
            &mut mat_c,
            ldc,
        );

        Matrix::from_data(c_nrows as usize, c_ncols as usize, mat_c)
    }

    // See the above comment about this string of identical functions with different names.  It's
    // for profiling but eventually you will need to get rid of them.

    pub unsafe fn mmul2_apple(
        &self,
        other: &Matrix,
        transpose_self: u8,
        transpose_other: u8,
    ) -> Result<Matrix, String> {
        // TODO check if the matrices can be multiplied
        let transpose_a = transpose_self;
        let transpose_b = transpose_other;
        let mat_a = &self.data;
        let mat_b = &other.data;
        let alpha = 1.;
        let beta = 1.;

        // m is nrows of Op(A) and nrows matrix C
        let m = if transpose_self == TRANSPOSE {
            // We're transposing self aka the a matrix
            self.ncols as i32
        } else {
            self.nrows as i32
        };

        // n is ncols of Op(b) and ncols of matrix c
        let n = if transpose_other == TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        // k is ncols of op(A), and nrows of op(B)
        let k = if transpose_self == TRANSPOSE && transpose_other == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_self == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_other == TRANSPOSE {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        };

        let c_nrows = m;
        let c_ncols = n;

        let mut mat_c = vec![0.; (c_nrows * c_ncols) as usize];

        let lda = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        let ldb = if transpose_other == NO_TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        let ldc = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        // let lda = self.nrows as i32;
        // let ldb = other.nrows as i32;
        // let ldc = self.nrows as i32;

        // DGEMM: C = alpha * op(A) * op(B) + beta * C.  Note: C gets overwritten with this.

        dgemm(
            transpose_a,
            transpose_b,
            m,
            n,
            k,
            alpha,
            mat_a,
            lda,
            mat_b,
            ldb,
            beta,
            &mut mat_c,
            ldc,
        );

        Matrix::from_data(c_nrows as usize, c_ncols as usize, mat_c)
    }

    pub unsafe fn mmul2_bonky(
        &self,
        other: &Matrix,
        transpose_self: u8,
        transpose_other: u8,
    ) -> Result<Matrix, String> {
        // TODO check if the matrices can be multiplied
        let transpose_a = transpose_self;
        let transpose_b = transpose_other;
        let mat_a = &self.data;
        let mat_b = &other.data;
        let alpha = 1.;
        let beta = 1.;

        // m is nrows of Op(A) and nrows matrix C
        let m = if transpose_self == TRANSPOSE {
            // We're transposing self aka the a matrix
            self.ncols as i32
        } else {
            self.nrows as i32
        };

        // n is ncols of Op(b) and ncols of matrix c
        let n = if transpose_other == TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        // k is ncols of op(A), and nrows of op(B)
        let k = if transpose_self == TRANSPOSE && transpose_other == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_self == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_other == TRANSPOSE {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        };

        let c_nrows = m;
        let c_ncols = n;

        let mut mat_c = vec![0.; (c_nrows * c_ncols) as usize];

        let lda = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        let ldb = if transpose_other == NO_TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        let ldc = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        // let lda = self.nrows as i32;
        // let ldb = other.nrows as i32;
        // let ldc = self.nrows as i32;

        // DGEMM: C = alpha * op(A) * op(B) + beta * C.  Note: C gets overwritten with this.

        dgemm(
            transpose_a,
            transpose_b,
            m,
            n,
            k,
            alpha,
            mat_a,
            lda,
            mat_b,
            ldb,
            beta,
            &mut mat_c,
            ldc,
        );

        Matrix::from_data(c_nrows as usize, c_ncols as usize, mat_c)
    }

    pub unsafe fn mmul2_carlisle(
        &self,
        other: &Matrix,
        transpose_self: u8,
        transpose_other: u8,
    ) -> Result<Matrix, String> {
        // TODO check if the matrices can be multiplied
        let transpose_a = transpose_self;
        let transpose_b = transpose_other;
        let mat_a = &self.data;
        let mat_b = &other.data;
        let alpha = 1.;
        let beta = 1.;

        // m is nrows of Op(A) and nrows matrix C
        let m = if transpose_self == TRANSPOSE {
            // We're transposing self aka the a matrix
            self.ncols as i32
        } else {
            self.nrows as i32
        };

        // n is ncols of Op(b) and ncols of matrix c
        let n = if transpose_other == TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        // k is ncols of op(A), and nrows of op(B)
        let k = if transpose_self == TRANSPOSE && transpose_other == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_self == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_other == TRANSPOSE {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        };

        let c_nrows = m;
        let c_ncols = n;

        let mut mat_c = vec![0.; (c_nrows * c_ncols) as usize];

        let lda = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        let ldb = if transpose_other == NO_TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        let ldc = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        // let lda = self.nrows as i32;
        // let ldb = other.nrows as i32;
        // let ldc = self.nrows as i32;

        // DGEMM: C = alpha * op(A) * op(B) + beta * C.  Note: C gets overwritten with this.

        dgemm(
            transpose_a,
            transpose_b,
            m,
            n,
            k,
            alpha,
            mat_a,
            lda,
            mat_b,
            ldb,
            beta,
            &mut mat_c,
            ldc,
        );

        Matrix::from_data(c_nrows as usize, c_ncols as usize, mat_c)
    }

    pub unsafe fn mmul2_dawn(
        &self,
        other: &Matrix,
        transpose_self: u8,
        transpose_other: u8,
    ) -> Result<Matrix, String> {
        // TODO check if the matrices can be multiplied
        let transpose_a = transpose_self;
        let transpose_b = transpose_other;
        let mat_a = &self.data;
        let mat_b = &other.data;
        let alpha = 1.;
        let beta = 1.;

        // m is nrows of Op(A) and nrows matrix C
        let m = if transpose_self == TRANSPOSE {
            // We're transposing self aka the a matrix
            self.ncols as i32
        } else {
            self.nrows as i32
        };

        // n is ncols of Op(b) and ncols of matrix c
        let n = if transpose_other == TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        // k is ncols of op(A), and nrows of op(B)
        let k = if transpose_self == TRANSPOSE && transpose_other == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_self == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_other == TRANSPOSE {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        };

        let c_nrows = m;
        let c_ncols = n;

        let mut mat_c = vec![0.; (c_nrows * c_ncols) as usize];

        let lda = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        let ldb = if transpose_other == NO_TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        let ldc = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        // let lda = self.nrows as i32;
        // let ldb = other.nrows as i32;
        // let ldc = self.nrows as i32;

        // DGEMM: C = alpha * op(A) * op(B) + beta * C.  Note: C gets overwritten with this.

        dgemm(
            transpose_a,
            transpose_b,
            m,
            n,
            k,
            alpha,
            mat_a,
            lda,
            mat_b,
            ldb,
            beta,
            &mut mat_c,
            ldc,
        );

        Matrix::from_data(c_nrows as usize, c_ncols as usize, mat_c)
    }

    pub unsafe fn mmul2_edgy(
        &self,
        other: &Matrix,
        transpose_self: u8,
        transpose_other: u8,
    ) -> Result<Matrix, String> {
        // TODO check if the matrices can be multiplied
        let transpose_a = transpose_self;
        let transpose_b = transpose_other;
        let mat_a = &self.data;
        let mat_b = &other.data;
        let alpha = 1.;
        let beta = 1.;

        // m is nrows of Op(A) and nrows matrix C
        let m = if transpose_self == TRANSPOSE {
            // We're transposing self aka the a matrix
            self.ncols as i32
        } else {
            self.nrows as i32
        };

        // n is ncols of Op(b) and ncols of matrix c
        let n = if transpose_other == TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        // k is ncols of op(A), and nrows of op(B)
        let k = if transpose_self == TRANSPOSE && transpose_other == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_self == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_other == TRANSPOSE {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        };

        let c_nrows = m;
        let c_ncols = n;

        let mut mat_c = vec![0.; (c_nrows * c_ncols) as usize];

        let lda = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        let ldb = if transpose_other == NO_TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        let ldc = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        // let lda = self.nrows as i32;
        // let ldb = other.nrows as i32;
        // let ldc = self.nrows as i32;

        // DGEMM: C = alpha * op(A) * op(B) + beta * C.  Note: C gets overwritten with this.

        dgemm(
            transpose_a,
            transpose_b,
            m,
            n,
            k,
            alpha,
            mat_a,
            lda,
            mat_b,
            ldb,
            beta,
            &mut mat_c,
            ldc,
        );

        Matrix::from_data(c_nrows as usize, c_ncols as usize, mat_c)
    }

    pub unsafe fn mmul2_funzilla(
        &self,
        other: &Matrix,
        transpose_self: u8,
        transpose_other: u8,
    ) -> Result<Matrix, String> {
        // TODO check if the matrices can be multiplied
        let transpose_a = transpose_self;
        let transpose_b = transpose_other;
        let mat_a = &self.data;
        let mat_b = &other.data;
        let alpha = 1.;
        let beta = 1.;

        // m is nrows of Op(A) and nrows matrix C
        let m = if transpose_self == TRANSPOSE {
            // We're transposing self aka the a matrix
            self.ncols as i32
        } else {
            self.nrows as i32
        };

        // n is ncols of Op(b) and ncols of matrix c
        let n = if transpose_other == TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        // k is ncols of op(A), and nrows of op(B)
        let k = if transpose_self == TRANSPOSE && transpose_other == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_self == TRANSPOSE {
            let ncols_op_a = self.nrows; // cos it's transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else if transpose_other == TRANSPOSE {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.ncols; // cos it's transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        } else {
            let ncols_op_a = self.ncols; // cos it's NOT transposed
            let nrows_op_b = other.nrows; // cos it's NOT transposed
            assert_eq!(ncols_op_a, nrows_op_b);

            ncols_op_a as i32
        };

        let c_nrows = m;
        let c_ncols = n;

        let mut mat_c = vec![0.; (c_nrows * c_ncols) as usize];

        let lda = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        let ldb = if transpose_other == NO_TRANSPOSE {
            other.nrows as i32
        } else {
            other.ncols as i32
        };

        let ldc = if transpose_self == NO_TRANSPOSE {
            self.nrows as i32
        } else {
            self.ncols as i32
        };

        // let lda = self.nrows as i32;
        // let ldb = other.nrows as i32;
        // let ldc = self.nrows as i32;

        // DGEMM: C = alpha * op(A) * op(B) + beta * C.  Note: C gets overwritten with this.

        dgemm(
            transpose_a,
            transpose_b,
            m,
            n,
            k,
            alpha,
            mat_a,
            lda,
            mat_b,
            ldb,
            beta,
            &mut mat_c,
            ldc,
        );

        Matrix::from_data(c_nrows as usize, c_ncols as usize, mat_c)
    }

    // > m <- matrix(c(10, 1, 1, 10), 2)
    // > m
    // [,1] [,2]
    // [1,]   10    1
    // [2,]    1   10
    // > chol(m)
    // [,1]      [,2]
    // [1,] 3.162278 0.3162278
    // [2,] 0.000000 3.1464265
    // > t(chol(m)) %*% chol(m)
    // [,1] [,2]
    // [1,]   10    1
    // [2,]    1   10
    pub unsafe fn chol(&self) -> Result<Matrix, String> {
        // TODO also check symmetric
        if self.nrows == self.ncols {
            let uplo = b'L';
            let n = self.ncols as i32;
            let mut chol = self.data.clone();
            let lda = self.nrows as i32;
            let mut info = 0;

            // Note that on output, only the lower triangular part of chol is modified.
            dpotrf(uplo, n, &mut chol, lda, &mut info);

            match info.cmp(&0) {
                Ordering::Greater => {
                    let msg = format!("the leading minor order of {} is not positive definite, and the factorization could not be completed", info);
                    Err(msg)
                }
                Ordering::Less => {
                    let msg = format!("the {}-th argument had an illegal value", info);
                    Err(msg)
                }
                Ordering::Equal => {
                    let mut result = Matrix::from_data(self.nrows, self.nrows, chol).unwrap();

                    // Zero out the part not in lower triangular.
                    for j in 0..result.ncols {
                        for i in 0..result.nrows {
                            if i < j {
                                result.set(i, j, 0.);
                            }
                        }
                    }

                    Ok(result)
                }
            }
        } else {
            Err("self is not a square matrix".to_string())
        }
    }

    /// Compute all eigenvalues & eigenvectors
    pub unsafe fn eigen(&self) -> Result<EigenOutput, String> {
        // TODO dsyevr

        if self.nrows == self.ncols {
            // both eigenvectors and eigenvalues
            let jobz = b'V';
            // compute all of them
            let range = b'A';
            // stores lower triangular matrix.
            let uplo = b'L';
            // order of matrix self
            let n = self.nrows as i32;
            let lda = n;
            let ldz = n;

            // dim of this should be lda, *
            let mut a = self.data.clone();

            // This will hold the result of the workspace query.
            let mut work = vec![-1.];
            let lwork = -1; // first set to -1 to do the workspace query

            // Another workspace query...not sure if i can do them both at the same time!
            let mut iwork = vec![-1];
            let liwork = -1;

            // These are not referenced because range is b'A'.  They specify different ways of
            // searching for eigenvalues.
            let vl = 0.;
            let vu = 0.;
            let il = 0;
            let iu = 0;

            // Absolute error tolerance to which each eigenvalue/vector is required.  May be
            // different from this in certain situations....see documentation.
            let abstol = 1e-06;

            // holds total number of eigen values found
            let mut m = -1; // this should be n after running since range is 'A'
                            // Stores eigenvalues in ascending order
            let mut w = vec![-1.; n as usize];
            // TODO describe this
            // dim will end up being (ldz, m) but m should be 'n' because range is 'A'
            // if jobz is 'V' (which it is here) and info = 0 (ie success) then the columns of z
            // contain the eigen vectors.  the j-th column of z holds the eigen vector for w(j)
            // eigenvalue.
            let mut z = vec![-1.; (ldz * n) as usize];

            // size is 2 * max(1, m) but m isn't known yet...so it's max should be n.  The i-th
            // eigenvector is nonzero only in elements isuppz(2i-1) through isuppz(2i).
            let mut isuppz = vec![-1; (2 * n) as usize]; // probably could be ldz * n as well?

            // 0: success, -i: the i-th parameter was bad, i: an internal error has occurred.
            let mut info = 0;

            // Computes selected eigenvalues and eigenvectors of real symmetric matrix using
            // Relatively Robust Representations.

            // The first call is the workspace query.
            dsyevr(
                jobz,
                range,
                uplo,
                n,
                &mut a,
                lda,
                vl,
                vu,
                il,
                iu,
                abstol,
                &mut m,
                &mut w,
                &mut z,
                ldz,
                &mut isuppz,
                &mut work,
                lwork,
                &mut iwork,
                liwork,
                &mut info,
            );

            // Check the results of the workspace query.
            let lwork = match info.cmp(&0) {
                Ordering::Equal => {
                    // If work[0] is still -1, then the workspace query failed.
                    approx::assert_abs_diff_ne!(work[0], -1., epsilon = TOL);

                    work[0].ceil() as i32
                }
                // TODO put in the real error message
                Ordering::Less => panic!("bad1"),
                // TODO put in the real error message
                Ordering::Greater => panic!("bad2"),
            };

            let liwork = match info.cmp(&0) {
                Ordering::Equal => {
                    // If work[0] is still -1, then the workspace query failed.
                    assert_ne!(iwork[0], -1);

                    iwork[0]
                }
                // TODO put in the real error message
                Ordering::Less => panic!("bad3"),
                // TODO put in the real error message
                Ordering::Greater => panic!("bad4"),
            };

            let mut work = vec![-1.; lwork as usize];
            let mut iwork = vec![-1; liwork as usize];

            // TODO better error message.
            // if info != 0, then something went wrong.
            assert_eq!(info, 0);

            // Now run the actual computation.
            dsyevr(
                jobz,
                range,
                uplo,
                n,
                &mut a,
                lda,
                vl,
                vu,
                il,
                iu,
                abstol,
                &mut m,
                &mut w,
                &mut z,
                ldz,
                &mut isuppz,
                &mut work,
                lwork,
                &mut iwork,
                liwork,
                &mut info,
            );

            // m should == n because we wanted all eigenvalues and eigenvectors.
            assert_eq!(m, n);

            // If it's not 0, then the computation failed.
            // TODO put in a real error message.
            assert_eq!(info, 0);

            Ok(EigenOutput {
                values: w,
                vectors: Matrix::from_data(n as usize, n as usize, z).unwrap(),
            })
        } else {
            // TODO also check for symmetry
            Err("the input was not a square matrix".to_string())
        }
    }

    /// Solves Ax = b, for x, where A is self, `b` is another Matrix.
    ///
    /// `self` is m X n
    pub unsafe fn solve(&self, other: &Matrix) -> Result<Matrix, String> {
        let trans = b'N';
        let m = self.nrows as i32;
        let n = self.ncols as i32;
        // aka num right hand sides
        let nrhs = other.ncols as i32;
        let mut a = self.data.clone();
        let mut b = other.data.clone();
        let lda = self.nrows as i32;
        let ldb = other.nrows as i32;

        // This will hold the result of the workspace query.
        let mut tmp = vec![0.];
        let lwork = -1;
        let mut info = 0;

        // The first call is the workspace query.  It will return with the proper value of lwork.
        dgels(
            trans, m, n, nrhs, &mut a, lda, &mut b, ldb, &mut tmp, lwork, &mut info,
        );

        // Check the result.
        let return_value = match info.cmp(&0) {
            Ordering::Equal => {
                // It was a success!  Now tmp[0] contains the min value of lwork for optimum
                // performance.
                let lwork = tmp[0].round() as i32;

                let mut work = vec![0.; lwork as usize];

                // Now we actaully run dgels.
                dgels(
                    trans, m, n, nrhs, &mut a, lda, &mut b, ldb, &mut work, lwork, &mut info,
                );

                // Need to match on info once more to check the output.
                let rv = match info.cmp(&0) {
                    Ordering::Equal => {
                        // was good
                        // TODO don't forget when testing to test three things at least...m>n, m<n, m=n.

                        // get the first n rows of b, and that's the output.
                        // the output size will be A.ncols X B.ncols
                        //
                        // TODO this is a roundabout way to get what I want.
                        let mut result = Matrix::zeros(self.ncols, other.ncols);

                        // b go through it and transfer it into result.
                        let tmp_b = Matrix::from_data(other.nrows, other.ncols, b.clone()).unwrap();

                        for j in 0..result.ncols() {
                            for i in 0..result.nrows() {
                                result.set(i, j, tmp_b.get(i, j));
                            }
                        }

                        Ok(result)
                    }
                    Ordering::Less => {
                        // Fail!
                        let msg =
                            format!("from dgels: the {}-th parameter had an illegal value", info);

                        Err(msg)
                    }
                    Ordering::Greater => {
                        // Also fail!
                        let msg = format!("from dgels: the {}-th diagonal element of the triangular factor of A is zero, so
        that A does not have full rank; the least squares solution could not be
        computed", info);

                        Err(msg)
                    }
                };

                rv
            }
            Ordering::Less => {
                // Fail!
                let msg = format!("from dgels: the {}-th parameter had an illegal value", info);

                Err(msg)
            }
            Ordering::Greater => {
                // Also fail!
                let msg = format!("from dgels: the {}-th diagonal element of the triangular factor of A is zero, so
        that A does not have full rank; the least squares solution could not be
        computed", info);

                Err(msg)
            }
        };

        return_value
    }

    pub fn center(&self, by_col: bool) -> Matrix {
        let mut result = Matrix::zeros(self.nrows, self.ncols);

        if by_col {
            let colmeans = self.colmeans();

            for j in 0..self.ncols {
                let mean = colmeans[j];

                for i in 0..self.nrows {
                    result.set(i, j, self.get(i, j) - mean);
                }
            }

            result
        } else {
            let rowmeans = self.rowmeans();

            for i in 0..self.nrows {
                let mean = rowmeans[i];

                for j in 0..self.ncols {
                    result.set(i, j, self.get(i, j) - mean);
                }
            }

            result
        }
    }

    pub fn colmeans(&self) -> Vec<f64> {
        let mut result = Vec::new();

        for j in 0..self.ncols {
            let col = self.col(j);

            let mean = col.iter().sum::<f64>() / col.len() as f64;

            result.push(mean);
        }

        result
    }

    // TODO test this
    pub fn rowmeans(&self) -> Vec<f64> {
        let mut result = Vec::new();

        for i in 0..self.nrows {
            let row = self.row(i);

            let mean = row.iter().sum::<f64>() / row.len() as f64;

            result.push(mean);
        }

        result
    }

    pub fn covariance(&self) -> Matrix {
        let mut result = Matrix::zeros(self.ncols, self.ncols);

        for j in 0..result.ncols() {
            for i in j..result.nrows() {
                if i == j {
                    result.set(i, j, var(self.col(j)));
                } else if i > j {
                    let cov = cov(self.col(i), self.col(j));

                    result.set(i, j, cov);
                    result.set(j, i, cov);
                }
            }
        }

        result
    }

    /// Gets all elements ij where i == j.
    ///
    /// So it works on square or not square!
    pub fn diag(&self) -> Vec<f64> {
        let mut diagonal = Vec::new();

        for j in 0..self.ncols {
            for i in 0..self.nrows {
                if i == j {
                    diagonal.push(self.get(i, j));
                }
            }
        }

        diagonal
    }

    /// Consumes `v`.
    pub fn from_diag(v: &[f64]) -> Matrix {
        let mut ary = Matrix::zeros(v.len(), v.len());

        for j in 0..ary.ncols() {
            ary.set(j, j, v[j]);
        }

        ary
    }

    /// Transpose.
    ///
    /// Avoid at all costs!  It's slow!  If these are taking a lot of time, add a better algorithm.
    pub fn transpose(&self) -> Matrix {
        let mut result = Matrix::zeros(self.ncols, self.nrows);

        for j in 0..result.ncols {
            for i in 0..result.nrows {
                result.set(i, j, self.get(j, i));
            }
        }

        result
    }
}

/// Variance of vector....uses N - 1 for the denominator to give unbiased estimator for iid
/// observations (See R).
pub fn var(v: &[f64]) -> f64 {
    let mean = v.iter().sum::<f64>() / v.len() as f64;

    v.iter().map(|&x| ((x - mean) as f64).powi(2)).sum::<f64>() / (v.len() - 1) as f64
}

// TODO return Result
pub fn cov(v1: &[f64], v2: &[f64]) -> f64 {
    let len = v1.len();
    assert_eq!(len, v2.len());

    let v1mu = v1.iter().sum::<f64>() / v1.len() as f64;
    let v2mu = v2.iter().sum::<f64>() / v2.len() as f64;

    // let (vv1, vv2) = (&v1[..len], &v2[..len]);
    //
    // let mut sum = 0.;
    // for i in 0..len {
    //     sum += (vv1[i] - v1mu) * (vv2[i] - v2mu)
    // }
    //
    // sum / (len - 1) as f64

    v1.iter()
        .zip(v2)
        .map(|(&x, &y)| (x - v1mu) * (y - v2mu))
        .sum::<f64>()
        / (v1.len() - 1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx;
    use std::f64;

    fn make_test_matrix() -> Matrix {
        Matrix::from_data(
            3,
            4,
            vec![
                1., 2., 3., 10., 20., 30., 100., 200., 300., 1000., 2000., 3000.,
            ],
        )
        .unwrap()
    }

    fn make_small_matrix() -> Matrix {
        Matrix::from_data(2, 2, vec![1., 2., 10., 20.]).unwrap()
    }

    #[test]
    fn get_gets_elements() {
        let m = make_test_matrix();

        for j in 0..m.ncols() {
            let multiplier: f64 = 10.0_f64.powi(j as i32);

            for i in 0..m.nrows() {
                approx::assert_abs_diff_eq!(
                    m.get(i, j),
                    (1. + i as f64) * multiplier,
                    epsilon = TOL
                );
            }
        }
    }

    #[test]
    fn get_unchecked_gets_elements_without_bounds_checking() {
        let m = make_test_matrix();

        for j in 0..m.ncols() {
            let multiplier: f64 = 10.0_f64.powi(j as i32);

            for i in 0..m.nrows() {
                unsafe {
                    approx::assert_abs_diff_eq!(
                        m.get_unchecked(i, j),
                        (1. + i as f64) * multiplier,
                        epsilon = TOL
                    );
                }
            }
        }
    }

    #[test]
    fn set_sets_elements() {
        let mut m = make_test_matrix();

        for j in 0..m.ncols() {
            for i in 0..m.nrows() {
                m.set(i, j, 0.);
            }
        }

        let expected = Matrix::zeros(m.nrows(), m.ncols());

        approx::assert_abs_diff_eq!(m.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn set_unchecked_sets_elements_without_bounds_checking() {
        let mut m = make_test_matrix();

        for j in 0..m.ncols() {
            for i in 0..m.nrows() {
                unsafe {
                    m.set_unchecked(i, j, 0.);
                }
            }
        }

        let expected = Matrix::zeros(m.nrows(), m.ncols());

        approx::assert_abs_diff_eq!(m.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn row_sums_sums_rows() {
        let m = make_test_matrix();
        let expected = vec![1111., 2222., 3333.];
        let actual = m.row_sums();

        approx::assert_abs_diff_eq!(&actual[..], &expected[..], epsilon = TOL);
    }

    #[test]
    fn col_sums_sums_colss() {
        let m = make_test_matrix();
        let expected = vec![6., 60., 600., 6000.];
        let actual = m.col_sums();

        approx::assert_abs_diff_eq!(&actual[..], &expected[..], epsilon = TOL);
    }

    #[test]
    fn transpose_transposes_a_matrix() {
        let a = Matrix::from_data(2, 3, vec![1., 2., 3., 4., 5., 6.]).unwrap();
        let actual = a.transpose();
        let expected = Matrix::from_data(3, 2, vec![1., 3., 5., 2., 4., 6.]).unwrap();

        approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn col_returns_a_column_slice() {
        let a = make_test_matrix();

        approx::assert_abs_diff_eq!(a.col(0), &vec![1., 2., 3.][..], epsilon = TOL);
        approx::assert_abs_diff_eq!(a.col(1), &vec![10., 20., 30.][..], epsilon = TOL);
        approx::assert_abs_diff_eq!(a.col(2), &vec![100., 200., 300.][..], epsilon = TOL);
        approx::assert_abs_diff_eq!(a.col(3), &vec![1000., 2000., 3000.][..], epsilon = TOL);
    }

    #[test]
    fn multiple_cols_returns_mutliple_columns() {
        let a = make_test_matrix();

        // Requesting a single row matches the .col() function.
        for j in 0..a.ncols() {
            approx::assert_abs_diff_eq!(a.multiple_cols(j, 1), a.col(j), epsilon = TOL);
        }

        approx::assert_abs_diff_eq!(
            a.multiple_cols(0, 2),
            &vec![1., 2., 3., 10., 20., 30.,][..],
            epsilon = TOL
        );

        approx::assert_abs_diff_eq!(
            a.multiple_cols(1, 3),
            &vec![10., 20., 30., 100., 200., 300., 1000., 2000., 3000.,][..],
            epsilon = TOL
        );

        approx::assert_abs_diff_eq!(a.multiple_cols(0, a.ncols()), a.data(), epsilon = TOL);
    }

    #[test]
    fn solve() {
        // This is the example from example_DGELS_colmajor.c of the LAPACK package.
        let a = Matrix::from_data(
            5,
            3,
            vec![1., 2., 3., 4., 5., 1., 3., 5., 2., 4., 1., 4., 2., 5., 3.],
        )
        .unwrap();
        let b = Matrix::from_data(
            5,
            2,
            vec![-10., 12., 14., 16., 18., -3., 14., 12., 16., 16.],
        )
        .unwrap();

        let expected = Matrix::from_data(3, 2, vec![2., 1., 1., 1., 1., 2.]).unwrap();

        let actual = unsafe { a.solve(&b).unwrap() };

        approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn eigen_gets_eigenvalues_and_vectors() {
        // This is covariance matrix of first 3 variables in iris data set.
        let a = Matrix::from_data(
            3,
            3,
            vec![0.69, -0.04, 1.27, -0.04, 0.19, -0.33, 1.27, -0.33, 3.12],
        )
        .unwrap();

        let actual = unsafe { a.eigen().unwrap() };
        let expected = EigenOutput {
            values: vec![0.06141744, 0.24687849, 3.69170407],
            vectors: Matrix::from_data(
                3,
                3,
                vec![
                    // First
                    -0.6499097,
                    0.6781939,
                    0.3430312,
                    // Middle
                    0.6528732,
                    0.7292474,
                    -0.2048286,
                    // Last
                    0.38906812,
                    -0.09083583,
                    0.91671962,
                ],
            )
            .unwrap(),
        };

        approx::assert_abs_diff_eq!(&actual.values[..], &expected.values[..], epsilon = TOL);
        approx::assert_abs_diff_eq!(
            actual.vectors.data(),
            expected.vectors.data(),
            epsilon = TOL
        );
    }

    #[test]
    fn mmul2_multiplies_two_matrices() {
        let a = Matrix::from_data(2, 2, vec![1., 2., 3., 4.]).unwrap();

        let a_x_a = unsafe { a.mmul2(&a, NO_TRANSPOSE, NO_TRANSPOSE).unwrap() };
        let at_x_a = unsafe { a.mmul2(&a, TRANSPOSE, NO_TRANSPOSE).unwrap() };
        let a_x_at = unsafe { a.mmul2(&a, NO_TRANSPOSE, TRANSPOSE).unwrap() };
        let at_x_at = unsafe { a.mmul2(&a, TRANSPOSE, TRANSPOSE).unwrap() };

        // Answers come from R
        approx::assert_abs_diff_eq!(a_x_a.data(), &vec![7., 10., 15., 22.][..], epsilon = TOL);
        approx::assert_abs_diff_eq!(at_x_a.data(), &vec![5., 11., 11., 25.][..], epsilon = TOL);
        approx::assert_abs_diff_eq!(a_x_at.data(), &vec![10., 14., 14., 20.][..], epsilon = TOL);
        approx::assert_abs_diff_eq!(at_x_at.data(), &vec![7., 15., 10., 22.][..], epsilon = TOL);
    }

    #[test]
    fn center_centers_cols_by_mean() {
        let a = make_test_matrix();

        let expected = Matrix::from_data(
            a.nrows(),
            a.ncols(),
            vec![
                -1., 0., 1., -10., 0., 10., -100., 0., 100., -1000., 0., 1000.,
            ],
        )
        .unwrap();

        let actual = a.center(true);

        approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn center_centers_rows_by_mean() {
        let a = Matrix::from_data(2, 3, vec![-1., -10., 0., 0., 1., 10.]).unwrap();

        let actual = a.center(false);

        approx::assert_abs_diff_eq!(actual.data(), a.data(), epsilon = TOL);
    }

    #[test]
    fn colmeans_gets_column_means() {
        let a = make_test_matrix();

        let expected = &vec![2., 20., 200., 2000.][..];

        let actual = &a.colmeans()[..];

        approx::assert_abs_diff_eq!(actual, expected, epsilon = TOL);
    }

    #[test]
    fn variance_gets_variance_of_slice() {
        let v = vec![1., 4., 3., 2., 7.];
        // R code:  var(c(1, 4, 3, 2, 7))
        let expected = 5.3;
        let actual = var(&v);

        approx::assert_abs_diff_eq!(actual, expected, epsilon = TOL);
    }

    #[test]
    fn covariance_gets_covariance_of_two_slices() {
        let x = vec![1., 4., 3., 2., 7.];
        let y = vec![10., 40., 30., 20., 70.];

        // R code:  var(c(1, 4, 3, 2, 7))
        let expected = 53.;
        let actual = cov(&x, &y);

        approx::assert_abs_diff_eq!(actual, expected, epsilon = TOL);
    }

    #[test]
    fn covariance_gets_cov_matrix() {
        let m = Matrix::zeros(2, 2);

        let actual = m.covariance();

        let expected = Matrix::from_data(2, 2, vec![0., 0., 0., 0.]).unwrap();

        approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);

        // another test

        let m = Matrix::from_data(2, 3, vec![1., 4., 3., 2., 7., 60.]).unwrap();

        let actual = m.covariance();

        let expected = Matrix::from_data(
            3,
            3,
            vec![4.5, -1.5, 79.5, -1.5, 0.5, -26.5, 79.5, -26.5, 1404.5],
        )
        .unwrap();

        approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn diag_gets_elements_on_the_diagonal() {
        let a = make_test_matrix();
        let expected = &vec![1., 20., 300.][..];
        let actual = a.diag();

        approx::assert_abs_diff_eq!(&actual[..], expected, epsilon = TOL);
    }

    #[test]
    fn to_diag_makes_a_vec_into_diagonal_array() {
        let v = vec![1., 2., 3.];
        let actual = Matrix::from_diag(&v);
        let expected = Matrix::from_data(3, 3, vec![1., 0., 0., 0., 2., 0., 0., 0., 3.]).unwrap();

        approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
    }

    #[test]
    fn row_gets_a_single_row_from_array() {
        let a = make_test_matrix();

        approx::assert_abs_diff_eq!(
            &a.row(0)[..],
            &vec![1., 10., 100., 1000.][..],
            epsilon = TOL
        );
        approx::assert_abs_diff_eq!(
            &a.row(1)[..],
            &vec![2., 20., 200., 2000.][..],
            epsilon = TOL
        );
        approx::assert_abs_diff_eq!(
            &a.row(2)[..],
            &vec![3., 30., 300., 3000.][..],
            epsilon = TOL
        );
    }

    mod equality {
        use super::*;

        #[test]
        fn eq_checks_equality() {
            let a = Matrix::from_data(2, 2, vec![1., 2., 3., 4.]).unwrap();
            let mut b = Matrix::from_data(2, 2, vec![1., 2., 3., 4.]).unwrap();

            assert!(a.eq(&b));

            b.set(0, 0, 100.);

            assert!(!a.eq(&b));
        }

        #[test]
        fn approx_eq_checks_equality_within_tolerance() {
            let a = Matrix::from_data(2, 2, vec![1., 2., 3., 4.]).unwrap();
            let mut b = Matrix::from_data(2, 2, vec![1.5, 2., 3., 4.]).unwrap();

            let tol = 0.5;

            assert!(a.approx_eq(&b, tol));

            // Bump the first element just beyond the tolerance.
            b.set(0, 0, 1.51);

            assert!(!a.approx_eq(&b, tol));
        }
    }

    // Testing matrix summary functions.
    mod summary {
        use super::*;

        #[test]
        fn sum_sums_all_elems() {
            let m = make_small_matrix();

            approx::assert_abs_diff_eq!(m.sum(), 33., epsilon = TOL);
        }
    }

    // Testing element wise operations on a single operand.
    mod element_wise_ops1 {
        use super::*;

        // For testing two operand element wise matrix operations.
        macro_rules! test_ewise_op1 {
            ($matrix_op:expr, $op:expr) => {{
                let a = make_small_matrix();

                let expected = Matrix::from_data(
                    a.nrows(),
                    a.ncols(),
                    vec![$op(1.), $op(2.), $op(10.), $op(20.)],
                )
                .unwrap();

                let actual = $matrix_op(&a);

                assert_eq!(actual.dim(), a.dim());
                approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
            }};
        }

        #[test]
        fn exp_does_element_wise_exp() {
            test_ewise_op1!(Matrix::ewise_exp, f64::exp);
        }

        #[test]
        fn ln_does_element_wise_natural_log() {
            test_ewise_op1!(Matrix::ewise_ln, f64::ln);
        }
    }

    // Testing element wise operations on two operands.
    mod element_wise_ops2 {
        use super::*;

        // For testing two operand element wise matrix operations.
        macro_rules! test_ewise_op2 {
            ($matrix_op:expr, $op:expr) => {{
                let a = make_small_matrix();
                let b = make_small_matrix();

                let expected = Matrix::from_data(
                    a.nrows(),
                    a.ncols(),
                    vec![$op(1., 1.), $op(2., 2.), $op(10., 10.), $op(20., 20.)],
                )
                .unwrap();

                let actual = $matrix_op(&a, &b);

                assert_eq!(actual.dim(), a.dim());
                approx::assert_abs_diff_eq!(actual.data(), expected.data(), epsilon = TOL);
            }};
        }

        #[test]
        fn eadd_does_elementwise_addition() {
            test_ewise_op2!(Matrix::ewise_add, Add::add);
        }

        #[test]
        fn ediv_does_elementwise_division() {
            test_ewise_op2!(Matrix::ewise_div, Div::div);
        }

        #[test]
        fn emul_does_elementwise_multiplication() {
            test_ewise_op2!(Matrix::ewise_mul, Mul::mul);
        }

        #[test]
        fn esub_does_elementwise_subtraction() {
            test_ewise_op2!(Matrix::ewise_sub, Sub::sub);
        }
    }
}
