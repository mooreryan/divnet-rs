/// Functions & macros for working with Vec<64> and slices elementwise.
use std::f64;
use std::ops::{Add, Div, Mul, Sub};

pub fn map(f: fn(f64) -> f64, xs: &[f64]) -> Vec<f64> {
    xs.iter().map(|&x| f(x)).collect()
}

// TODO: switch to returning Result rather than panicking
pub fn zipmap(f: fn(f64, f64) -> f64, xs: &[f64], ys: &[f64]) -> Vec<f64> {
    if xs.len() != ys.len() {
        panic!("xs and ys must be the same length!");
    }

    xs.iter().zip(ys).map(|(&x, &y)| f(x, y)).collect()
}

pub fn sum(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>()
}

macro_rules! define_unary_ops {
    ($($mod:ident::$op:ident),*) => {$(
        pub fn $op(xs: &[f64]) -> Vec<f64> {
            map($mod::$op, xs)
        }
    )*}
}

macro_rules! define_binary_ops {
    ($($mod:ident::$op:ident),*) => {$(
        pub fn $op(xs: &[f64], ys: &[f64]) -> Vec<f64> {
            zipmap($mod::$op, xs, ys)
        }
    )*}
}

define_unary_ops!(f64::ln, f64::exp);

define_binary_ops!(Add::add, Div::div, Mul::mul, Sub::sub);

#[cfg(test)]
mod tests {
    use super::*;
    use approx;
    use std::f64;

    const TOL: f64 = 1e-5;

    macro_rules! test_unary_ops {
        ($($test_name:ident => $op:ident => $mod:ident::$orig_op:ident),*) => {$(
            #[test]
            fn $test_name() {
                let a = vec![1., 2.];

                let expected = vec![$mod::$op(a[0]), $mod::$op(a[1])];

                let actual = super::$op(&a);

                approx::assert_abs_diff_eq!(&actual[..], &expected[..], epsilon = TOL);
            }
        )*}
    }

    macro_rules! test_binary_ops {
        ($($test_name:ident => $op:ident => $trait:ident::$orig_op:ident),*) => {$(
            #[test]
            fn $test_name() {
                let a = vec![1., 2.];
                let b = vec![10., 20.];

                let expected = vec![$trait::$orig_op(a[0], b[0]), $trait::$orig_op(a[1], b[1])];

                let actual = super::$op(&a, &b);

                approx::assert_abs_diff_eq!(&actual[..], &expected[..], epsilon = TOL);
            }
        )*}
    }

    test_unary_ops!(test_ln => ln => f64::ln);

    test_binary_ops!(test_add => add => Add::add, test_div => div => Div::div, test_mul => mul => Mul::mul, test_sub => sub => Sub::sub);

    #[test]
    fn map_maps_a_fn_onto_a_vec() {
        let xs = vec![1., 2.];

        let actual = map(f64::ln, &xs);
        let expected = vec![f64::ln(1.), f64::ln(2.)];

        approx::assert_abs_diff_eq!(&actual[..], &expected[..], epsilon = TOL);
    }

    #[test]
    fn zipmap_maps_a_fn_onto_two_vecs_lockstep() {
        let xs = vec![1., 2.];
        let ys = vec![10., 20.];

        let actual = zipmap(Add::add, &xs, &ys);
        let expected = vec![11., 22.];

        approx::assert_abs_diff_eq!(&actual[..], &expected[..], epsilon = TOL);
    }

    #[test]
    #[should_panic]
    fn zipmap_panics_if_vecs_are_different_lengths() {
        let xs = vec![1., 2.];
        let ys = vec![1., 2., 3.];

        zipmap(Add::add, &xs, &ys);
    }

    #[test]
    fn sum_sums_elements() {
        let xs = vec![1., 2., 3.];

        let actual = sum(&xs);
        let expected = 6.;

        approx::assert_abs_diff_eq!(actual, expected, epsilon = TOL);
    }
}
