use blstrs::Bls12;
use blstrs::G1Affine;
use blstrs::G1Projective;
use blstrs::G2Affine;
use blstrs::G2Projective;
use blstrs::Gt;
use blstrs::Scalar;

use group::ff::Field;
use group::prime::PrimeCurveAffine;
use group::Curve;
use group::Group;
use pairing::Engine;

use nalgebra::DMatrix;

use crate::bls_elements::*;

pub fn gen_Γ_linear_fc(
    cm: G1Projective,
    ck: (Vec<G1Projective>, Vec<G2Projective>),
) -> DMatrix<BlsElement> {
    let cm_b = BlsElement::G1Affine(cm.to_affine());
    // let g1_1 = BlsElement::G1Affine(G1Affine::identity());
    // let g2_1 = BlsElement::G2Affine(G2Affine::identity());
    let u2 = BlsElement::G2Affine(ck.1[0].to_affine());

    let Γ = DMatrix::from_row_slice(1, 2, vec![cm_b, u2].as_slice());

    Γ
}

pub fn mul_Γ(A: DMatrix<BlsElement>, Γ: DMatrix<BlsElement>) -> DMatrix<BlsElement> {
    let width = Γ.ncols();
    let mut res: Vec<BlsElement> = vec![];
    for i in 0..width {
        let mut element = A[0] * Γ[i * A.len()];
        for j in 1..A.len() {
            element = element + A[j] * Γ[i * A.len() + j];
        }
        res.push(element);
    }
    DMatrix::from_row_slice(width, 1, &res.as_slice())
}

pub fn mul_θ(A: DMatrix<BlsElement>, θ: DMatrix<BlsElement>) -> BlsElement {
    let mut element = A[0] * θ[0];
    for i in 1..A.len() {
        element = element + A[i] * θ[i];
    }
    element
}

pub fn mul_λ(A: DMatrix<BlsElement>, λ: DMatrix<BlsElement>) -> BlsElement {
    let mut element = A[0] * λ[0];
    for i in 1..A.len() {
        element = element + A[i] * λ[i];
    }
    element
}

pub fn gen_hash_keys(Γ: DMatrix<BlsElement>) -> (DMatrix<BlsElement>, DMatrix<BlsElement>) {
    let mut rng = rand::thread_rng();
    let a_n = DMatrix::from_fn(Γ.nrows(), 1, |_r, _c| {
        BlsElement::Scalar(Scalar::random(&mut rng))
    });
    let hp = mul_Γ(a_n.clone(), Γ);

    (a_n, hp)
}

pub fn gen_θ_linear_fc(
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    y: Scalar,
) -> DMatrix<BlsElement> {
    let θ = Bls12::pairing(
        &ck.0[0].to_affine(),
        &(ck.1[ck.1.len() - 1] * y).to_affine(),
    );
    DMatrix::from_vec(1, 1, vec![BlsElement::Gt(θ)])
}

pub fn gen_verifier_hash(hash_key: DMatrix<BlsElement>, θ: DMatrix<BlsElement>) -> Gt {
    let res = mul_θ(hash_key, θ);
    if let BlsElement::Gt(result) = res {
        return result;
    }
    panic!("We did not get back a Gt element");
}

pub fn gen_λ_linear_fc(element1: G2Projective, element2: G1Projective) -> DMatrix<BlsElement> {
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
