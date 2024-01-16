use bls12_381::Bls12;
use bls12_381::G1Affine;
use bls12_381::G1Projective;
use bls12_381::G2Affine;
use bls12_381::G2Projective;
use bls12_381::Gt;
use bls12_381::Scalar;

use ff::Field;
use group::Curve;
use pairing::Engine;

use nalgebra::DMatrix;

use crate::bls_elements::*;
use crate::linear_fc::*;

pub fn gen_gamma_linear_fc(cm: &G1Projective, ck: &CommitmentKey) -> DMatrix<BlsElement> {
    let cm_b = BlsElement::G1Affine(cm.to_affine());
    let u2 = BlsElement::G2Affine(ck.u2[0].to_affine());

    let gamma = DMatrix::from_row_slice(1, 2, vec![cm_b, u2].as_slice());

    gamma
}

pub fn mul_gamma(A: DMatrix<BlsElement>, gamma: DMatrix<BlsElement>) -> DMatrix<BlsElement> {
    let width = gamma.ncols();
    let mut res: Vec<BlsElement> = vec![];
    for i in 0..width {
        let mut element = A[0] * gamma[i * A.len()];
        for j in 1..A.len() {
            element = element + A[j] * gamma[i * A.len() + j];
        }
        res.push(element);
    }
    DMatrix::from_row_slice(width, 1, &res.as_slice())
}

pub fn mul_theta(A: DMatrix<BlsElement>, theta: DMatrix<BlsElement>) -> BlsElement {
    let mut element = A[0] * theta[0];
    for i in 1..A.len() {
        element = element + A[i] * theta[i];
    }
    element
}

pub fn mul_lambda(A: DMatrix<BlsElement>, lambda: DMatrix<BlsElement>) -> BlsElement {
    let mut element = A[0] * lambda[0];
    for i in 1..A.len() {
        element = element + A[i] * lambda[i];
    }
    element
}

pub fn gen_hash_keys(gamma: DMatrix<BlsElement>) -> (DMatrix<BlsElement>, DMatrix<BlsElement>) {
    let mut rng = rand::thread_rng();
    let a_n = DMatrix::from_fn(gamma.nrows(), 1, |_r, _c| {
        BlsElement::Scalar(Scalar::random(&mut rng))
    });
    let hp = mul_gamma(a_n.clone(), gamma);

    (a_n, hp)
}

pub fn gen_theta_linear_fc(ck: &CommitmentKey, y: Scalar) -> DMatrix<BlsElement> {
    let theta = Bls12::pairing(
        &ck.u1[0].to_affine(),
        &(ck.u2[ck.u2.len() - 1] * y).to_affine(),
    );
    DMatrix::from_vec(1, 1, vec![BlsElement::Gt(theta)])
}

pub fn gen_verifier_hash(hash_key: DMatrix<BlsElement>, theta: DMatrix<BlsElement>) -> Gt {
    let res = mul_theta(hash_key, theta);
    if let BlsElement::Gt(result) = res {
        return result;
    }
    panic!("We did not get back a Gt element");
}

pub fn gen_lambda_linear_fc(element1: G2Projective, element2: G1Projective) -> DMatrix<BlsElement> {
    DMatrix::from_row_slice(
        2,
        1,
        vec![
            BlsElement::G2Affine(element1.to_affine()),
            BlsElement::G1Affine(element2.to_affine()),
        ]
        .as_slice(),
    )
}
