pub mod bls_elements;
pub mod linear_fc;

use bitvec::prelude::Msb0;
use bitvec::view::BitViewSized;
use blstrs::{G1Projective, G2Projective, Gt, Scalar};
use group::{ff::Field, Group};
use linear_fc::{FunctionCoeffs, FunctionVars};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

pub fn encrypt_bit(
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    cm: G1Projective,
    beta: FunctionCoeffs,
    y: Scalar,
    message: bool,
) -> Result<(bool, Scalar, Gt, Vec<u8>), &'static str> {
    let mut rng = rand::thread_rng();

    let hash_key = vec![Scalar::random(&mut rng); 1];
    let projected_hash_key = hash_key[0];
    let m = linear_fc::righthandside(ck.clone(), y, beta.len() as u64);
    let projected_hash_key_2 = m * hash_key[0];
    let theta = linear_fc::lefthandside(cm, ck.clone(), beta.clone(), beta.len() as u64);
    let hash = theta * hash_key[0];

    let mut r: Vec<u8> = vec![];
    for _ in 0..linear_fc::L {
        let rand_byte = rand::random::<u8>();
        r.push(rand_byte);
    }
    let ciphertext = encrypt_message(message, hash, r.clone());

    Ok((ciphertext, projected_hash_key, projected_hash_key_2, r))
}

pub fn decrypt_bit(
    ciphertext: bool,
    u: Scalar,
    beta: FunctionCoeffs,
    projected_hash_key: Scalar,
    projected_hash_key_2: Gt,
    r: Vec<u8>,
    d: (FunctionVars, Scalar),
) -> bool {
    let dl_of_opening = linear_fc::dl_of_opening(u, d.clone(), beta.clone());
    let projected_hash =
        projected_hash_key_2 + (Gt::generator() * (dl_of_opening * projected_hash_key));
    encrypt_message(ciphertext, projected_hash, r)
}

fn encrypt_message(message: bool, hashkey: Gt, rand_bytes: Vec<u8>) -> bool {
    let mut hasher = Shake256::default();
    hasher.update(hashkey.to_string().into_bytes().as_slice());
    let mut reader = hasher.finalize_xof();
    let mut hashedkey = [0u8; linear_fc::L as usize];
    reader.read(&mut hashedkey);

    let mut res_bit = false;
    for i in 0..rand_bytes.len() {
        let rand_bits: bitvec::array::BitArray<u8, Msb0> = rand_bytes[i].into_bitarray();
        let h_bits: bitvec::array::BitArray<u8, Msb0> = hashedkey[i].into_bitarray();
        for j in 0..8 {
            res_bit = res_bit ^ (rand_bits[j] ^ h_bits[j]);
        }
    }
    let result = res_bit ^ message;
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encryption_decryption(message: bool) -> bool {
        let witness_length = 2;
        let (ck, alpha, beta, y, u) = linear_fc::setup_random(witness_length);
        let (cm, d) = linear_fc::commit(ck.clone(), alpha);

        let (ciphertext, projected_hash_key, projected_hash_key_2, r) =
            encrypt_bit(ck.clone(), cm, beta.clone(), y, message).unwrap();
        decrypt_bit(
            ciphertext,
            u,
            beta,
            projected_hash_key,
            projected_hash_key_2,
            r,
            d,
        )
    }

    #[test]
    fn encryption_decryption_success() {
        for _ in 0..128 {
            // we're only encrypting a bit, which can be easily randomly correct, so to properly test we should run this often
            let message = rand::random::<bool>();
            assert_eq!(encryption_decryption(message), message);
        }
    }
}
