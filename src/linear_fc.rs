use ark_bls12_381::Fr as ScalarField;
use ark_bls12_381::{Bls12_381, Fq12, G1Projective, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::Group;
use ark_ff::Field;
use ark_std::{UniformRand, Zero};
use rand::thread_rng;
use std::iter;

#[derive(Debug, Clone)]
pub struct CommitmentKey {
    pub u1: Vec<G1Projective>,
    pub u2: Vec<G2Projective>,
}

pub fn setup_unsafe(n: u64) -> CommitmentKey {
    let mut rng = rand::thread_rng();
    let u = ScalarField::rand(&mut rng);

    let mut u1: Vec<G1Projective> = (0..2 * n)
        .map(|j| G1Projective::generator() * u.pow(&[j + 1, 0, 0, 0]))
        .collect();
    // set trapdoor value to zero so opening proof works correctly
    u1[n as usize] = G1Projective::zero();

    let u2: Vec<G2Projective> = (0..n)
        .map(|j| G2Projective::generator() * u.pow(&[j + 1, 0, 0, 0]))
        .collect();

    CommitmentKey { u1, u2 }
}

pub fn compute_func(x: &Vec<ScalarField>, beta: &Vec<ScalarField>) -> ScalarField {
    x.iter().zip(beta).map(|(a, b)| a * b).sum()
}

pub fn commit(ckey: &CommitmentKey, x: &Vec<ScalarField>) -> (G1Projective, ScalarField) {
    let mut rng = rand::thread_rng();
    let r = ScalarField::rand(&mut rng);
    let sum: G1Projective = iter::once(G1Projective::generator() * r)
        .chain(x.iter().zip(ckey.u1.iter()).map(|(&a, &u)| u * a))
        .sum();
    (sum, r)
}

pub fn open(
    ckey: &CommitmentKey,
    x: &Vec<ScalarField>,
    r: ScalarField,
    beta: &Vec<ScalarField>,
) -> G1Projective {
    assert_eq!(x.len(), beta.len());

    let n = x.len();
    let mut opening = G1Projective::zero();
    for i in 0..n {
        // assumes that u1[n] = 0
        let wi: G1Projective = ckey.u1[n - i - 1] * r
            + x.iter()
                .zip(ckey.u1.iter().skip(n - i))
                .map(|(&a, &u)| u * a)
                .sum::<G1Projective>();
        opening += wi * beta[i];
    }
    opening
}

pub fn verify(
    ckey: &CommitmentKey,
    cm: &G1Projective,
    op: &G1Projective,
    beta: &Vec<ScalarField>,
    y: ScalarField,
) -> bool {
    let n = beta.len();
    let linear_comb: G2Projective = beta
        .iter()
        .zip(ckey.u2.iter().rev())
        .map(|(&b, &u)| u * b)
        .sum();

    let pairing1 = Bls12_381::pairing(cm, linear_comb);
    let pairing2 = Bls12_381::pairing(op, G2Projective::generator());
    let pairing3 = Bls12_381::pairing(&ckey.u1[0], &(ckey.u2[n - 1])) * y;

    pairing1 == pairing2 + pairing3
}

#[cfg(test)]
mod test {

    use super::*;

    fn do_test_fc(should_succeed: bool) -> bool {
        let n = 2;
        let ckey = setup_unsafe(n);
        let x = vec![ScalarField::from(1), ScalarField::from(2)];
        let beta = vec![ScalarField::from(3), ScalarField::from(2)];
        let y = compute_func(&x, &beta);
        let (commit, r) = commit(&ckey, &x);
        let mut op = open(&ckey, &x, r, &beta);

        if !should_succeed {
            let mut rng = ark_std::test_rng();
            op = op * ScalarField::rand(&mut rng);
        }

        verify(&ckey, &commit, &op, &beta, y)
    }

    #[test]
    fn fc_success() {
        assert_eq!(do_test_fc(true), true);
    }

    #[test]
    fn fc_failure() {
        assert_eq!(do_test_fc(false), false);
    }
}
