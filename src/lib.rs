mod bls_elements;
mod encrypt;
mod linear_fc;

use bls12_381::G1Projective;
use bls12_381::G2Projective;
use bls12_381::Gt;
use bls12_381::Scalar;
use bls_elements::BlsElement;

use bitvec::prelude::Msb0;
use bitvec::view::BitViewSized;
use nalgebra::DMatrix;
use sha3::digest::ExtendableOutput;
use sha3::digest::Update;
use sha3::digest::XofReader;
use sha3::Shake256;

use ff::Field;
use group::Curve;

pub fn encrypt_bit(
    message: bool,
    ck: linear_fc::CommitmentKey,
    x: Vec<Scalar>,
    y: Scalar,
) -> Result<(bool, DMatrix<BlsElement>, Scalar, Vec<u8>), &'static str> {
    let (cm, r) = linear_fc::commit(&ck, &x);

    let gamma = encrypt::gen_gamma_linear_fc(&cm, &ck);
    let (hash_key, projected_hash_key) = encrypt::gen_hash_keys(gamma);
    let theta = encrypt::gen_theta_linear_fc(&ck, y);
    let hash = encrypt::gen_verifier_hash(hash_key, theta);
    let mut r2: Vec<u8> = vec![];
    for _ in 0..linear_fc::LENGTH {
        let rand_byte = rand::random::<u8>();
        r2.push(rand_byte);
    }
    let ciphertext = encrypt_message(message, hash, r2.clone());

    Ok((ciphertext, projected_hash_key, r, r2))
}

pub fn decrypt_bit(
    ciphertext: bool,
    ck: linear_fc::CommitmentKey,
    x: Vec<Scalar>,
    beta: Vec<Scalar>,
    r: Scalar,
    projected_hash_key: DMatrix<BlsElement>,
    r2: Vec<u8>,
) -> bool {
    let opening = linear_fc::open(&ck, &x, r, &beta);

    // diverging from the WE_FC paper's implementation
    let ADJUSTMENT_3 = ck.u1[0].to_affine() * beta[1] * x[0]; // G^(u*alpha1*beta2)

    let lambda = encrypt::gen_lambda_linear_fc(
        ck.u2[1] * beta[0] + ck.u2[0] * beta[1],
        opening - ADJUSTMENT_3,
    );
    let res = projected_hash_key[0] * lambda[0] - projected_hash_key[1] * lambda[1];
    if let BlsElement::Gt(projected_hash) = res {
        let plaintext = encrypt_message(ciphertext, projected_hash, r2.clone());
        plaintext
    } else {
        panic!("mul_lambda did not return an element from Gt");
    }
}

fn encrypt_message(message: bool, hashkey: Gt, rand_bytes: Vec<u8>) -> bool {
    let mut hasher = Shake256::default();
    hasher.update(hashkey.to_string().into_bytes().as_slice());
    let mut reader = hasher.finalize_xof();
    let mut hashedkey = [0u8; linear_fc::LENGTH as usize];
    reader.read(&mut hashedkey);

    let mut res_bit = false;
    for i in 0..rand_bytes.len() {
        let rand_bits: bitvec::array::BitArray<u8, Msb0> = rand_bytes[i].into_bitarray();
        let h_bits: bitvec::array::BitArray<u8, Msb0> = hashedkey[i].into_bitarray();
        for j in 0..8 {
            res_bit = res_bit ^ (rand_bits[j] & h_bits[j]);
        }
    }
    let result = res_bit ^ message;
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encryption_decryption(should_succeed: bool, message: bool) -> bool {
        let witness_length = 2;
        let mut ck = linear_fc::setup_unsafe(witness_length);
        let mut x = vec![Scalar::one(), Scalar::one().double()];
        let mut beta = vec![Scalar::one().double(), Scalar::one()];
        let mut y = linear_fc::compute_func(&x, &beta);
        // println!("x = {:?}", x);
        // println!("beta = {:?}", beta);
        // println!("y = {:?}", y);
        let (ciphertext, projected_hash_key, r, r2) =
            encrypt_bit(message, ck.clone(), x.clone(), y).unwrap();
        if !should_succeed {
            y = y * y;
        }
        decrypt_bit(ciphertext, ck, x, beta, r, projected_hash_key, r2)
    }

    #[test]
    fn encryption_decryption_success() {
        let message = true;
        for _ in 0..128 {
            // we're only encrypting a bit, which can be easily randomly correct, so to properly test we should run this often
            assert_eq!(encryption_decryption(true, message), message);
        }
        // assert_eq!(encryption_decryption(true, message), message);
    }
}
