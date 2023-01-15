mod bls_elements;
mod linear_fc;
mod sphf;

use bls_elements::BlsElement;
use blstrs::Gt;
use blstrs::G1Projective;
use blstrs::G2Projective;
use blstrs::Scalar;

use nalgebra::DMatrix;
use group::Curve;
use sha3::Shake256;
use sha3::digest::Update;
use sha3::digest::ExtendableOutput;
use sha3::digest::XofReader;
use bitvec::prelude::Msb0;
use bitvec::view::BitViewSized;

pub fn encrypt_bit(message: bool, ck: (Vec<G1Projective>, Vec<G2Projective>), α: Vec<Scalar>, y: Scalar) -> Result<(bool, DMatrix<BlsElement>, Vec<u8>, (Vec<Scalar>, Scalar)), &'static str> {
    let (cm, d) = linear_fc::commit(ck.clone(), α.clone());

    let Γ = sphf::gen_Γ_linear_fc(cm, ck.clone());
    let (hash_key, projected_hash_key) = sphf::gen_hash_keys(Γ);
    let θ = sphf::gen_θ_linear_fc(ck.clone(), y);
    let hash = sphf::gen_verifier_hash(hash_key, θ);
    let mut r: Vec<u8> = vec![];
    for _ in 0..linear_fc::L {
        let rand_byte = rand::random::<u8>();
        r.push(rand_byte);
    }
    let ciphertext = encrypt_message(message, hash, r.clone());

    Ok((ciphertext, projected_hash_key, r, d))
}

pub fn decrypt_bit(ciphertext: bool, ck: (Vec<G1Projective>, Vec<G2Projective>), α: Vec<Scalar>, β: Vec<Scalar>, projected_hash_key: DMatrix<BlsElement>, r: Vec<u8>, d: (Vec<Scalar>, Scalar)) -> bool {
    let opening = linear_fc::open(ck.clone(), d, β.clone());

    // diverging from the WE_FC paper's implementation
    let ADJUSTMENT_3 = ck.0[0].to_affine()*β[1]*α[0]; // G^(u*α1*β2)

    let λ = sphf::gen_λ_linear_fc(
        ck.1[1]*β[0] + ck.1[0]*β[1],
        opening - ADJUSTMENT_3,
    );
    let res = projected_hash_key[0]*λ[0] - projected_hash_key[1]*λ[1];
    if let BlsElement::Gt(projected_hash) = res {
        let plaintext = encrypt_message(ciphertext, projected_hash, r.clone());
        plaintext
    } else {
        panic!("mul_λ did not return an element from Gt");
    }
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
        let (ck, α, β, mut y) = linear_fc::setup_random(witness_length);
        let (ciphertext, projected_hash_key, r, d) = encrypt_bit(message, ck.clone(), α.clone(), y).unwrap();
        if !should_succeed {
            y = y*y;
        }
        decrypt_bit(ciphertext, ck, α, β, projected_hash_key, r, d)
    }

    #[test]
    fn encryption_decryption_success() {
        let message = true;
        for _ in 0..128 { // we're only encrypting a bit, which can be easily randomly correct, so to properly test we should run this often
            assert_eq!(encryption_decryption(true, message), message);            
        }
    }
}
