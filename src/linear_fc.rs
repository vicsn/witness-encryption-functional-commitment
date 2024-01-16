use blstrs::Bls12;
use blstrs::G1Projective;
use blstrs::G2Affine;
use blstrs::G2Projective;
use blstrs::Scalar;

use group::ff::Field;
use group::prime::PrimeCurveAffine;
use group::Curve;
use group::Group;
use pairing::Engine;
use std::iter;

pub const LENGTH: u64 = 256; // Length of hashkey and random bits used to encrypt the message

#[derive(Debug, Clone)]
pub struct CommitmentKey {
    pub u1: Vec<G1Projective>,
    pub u2: Vec<G2Projective>,
}

pub fn setup_unsafe(n: u64) -> CommitmentKey {
    let mut rng = rand::thread_rng();
    let u = Scalar::random(&mut rng);

    let mut u1: Vec<G1Projective> = (0..2 * n)
        .map(|j| G1Projective::generator() * u.pow_vartime(&[j + 1]))
        .collect();
    // set trapdoor value to zero so opening proof works correctly
    u1[n as usize] = G1Projective::identity();

    let u2: Vec<G2Projective> = (0..n)
        .map(|j| G2Projective::generator() * u.pow_vartime(&[j + 1]))
        .collect();

    CommitmentKey { u1, u2 }
}

pub fn compute_func(x: &Vec<Scalar>, beta: &Vec<Scalar>) -> Scalar {
    x.iter().zip(beta).map(|(a, b)| a * b).sum()
}

pub fn commit(ck: &CommitmentKey, x: &Vec<Scalar>) -> (G1Projective, Scalar) {
    let mut rng = rand::thread_rng();
    let r = Scalar::random(&mut rng);
    let sum: G1Projective = iter::once(G1Projective::generator() * r)
        .chain(x.iter().zip(ck.u1.iter()).map(|(&a, u)| u * a))
        .sum();
    (sum, r)
}

pub fn open(ck: &CommitmentKey, x: &Vec<Scalar>, r: Scalar, beta: &Vec<Scalar>) -> G1Projective {
    assert_eq!(x.len(), beta.len());

    let n = x.len();
    let mut opening = G1Projective::identity();
    for i in 0..n {
        // assumes that u1[n] = 0
        let wi: G1Projective = ck.u1[n - i - 1] * r
            + x.iter()
                .zip(ck.u1.iter().skip(n - i))
                .map(|(&a, &u)| u * a)
                .sum::<G1Projective>();
        opening += wi * beta[i];
    }
    opening
}

pub fn verification_pairings(
    ck: &CommitmentKey,
    cm: &G1Projective,
    op: &G1Projective,
    beta: &Vec<Scalar>,
    y: Scalar,
) -> bool {
    let n = beta.len();
    let linear_comb: G2Projective = beta
        .iter()
        .zip(ck.u2.iter().rev())
        .map(|(&b, &u)| u * b)
        .sum();

    let pairing1 = Bls12::pairing(&cm.to_affine(), &linear_comb.to_affine());
    let pairing2 = Bls12::pairing(&op.to_affine(), &G2Affine::generator());
    let pairing3 = Bls12::pairing(&ck.u1[0].to_affine(), &(ck.u2[n - 1].to_affine())) * y;

    pairing1 == pairing2 + pairing3
}

#[cfg(test)]
mod test {

    use super::*;

    fn do_test_fc(should_succeed: bool) -> bool {
        let n = 2;
        let ck = setup_unsafe(n);
        let x = vec![Scalar::one(), Scalar::one().double()];
        let beta = vec![Scalar::one().double(), Scalar::one()];
        let y = compute_func(&x, &beta);
        let (cm, r) = commit(&ck, &x);
        let mut op = open(&ck, &x, r, &beta);

        if !should_succeed {
            let rng = rand::thread_rng();
            op = op * Scalar::random(rng);
        }

        verification_pairings(&ck, &cm, &op, &beta, y)
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
