use rand::distributions::WeightedIndex;
use rand::prelude::*;

pub fn multinomial<R: Rng>(mut rng: &mut R, _n: usize, size: usize, probs: &[f64]) -> Vec<f64> {
    let prob_sum = probs.iter().sum::<f64>();

    approx::assert_abs_diff_eq!(prob_sum, 1., epsilon = 1e-6);

    let mut draws = vec![0.; probs.len()];

    let mut min = 2.;
    for &x in probs.iter() {
        assert!(x <= 1.);
        if x < min {
            min = x;
        }
    }
    assert!(min > 0.);
    assert!(min <= 1.);
    let weights: Vec<usize> = probs.iter().map(|&x| (x / min).round() as usize).collect();

    let dist = WeightedIndex::new(&weights).unwrap();

    for _ in 0..size {
        let i = dist.sample(&mut rng);

        // increment this taxa
        draws[i] += 1.;
    }

    draws
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx;
    use std::f64;

    #[test]
    fn test_multinomial() {
        let probs = vec![0.25, 0.25, 0.5];

        // Do a lot of draws to get the accuracy up.
        let size = 100_000 as usize;

        // Todo change this to a seeded rng to test the output.
        let mut rng = thread_rng();

        let draws = multinomial(&mut rng, 1, size, &probs);

        let actual_probs: Vec<f64> = draws.iter().map(|&count| count / size as f64).collect();

        approx::assert_abs_diff_eq!(&actual_probs[..], &probs[..], epsilon = 0.05);
    }
}
