use ark_bls12_381::Fr as ScalarField;
use ark_bls12_381::{Bls12_381, Fq12, G1Projective, G2Projective};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::Group;
use ark_ff::Field;
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::UniformRand;

use rand::{thread_rng, Rng};

use crate::linear_fc::*;

pub struct Ciphertext {
    pub proj_key: G2Projective,
    pub rand_bytes: Vec<u8>,
    pub ciphertext: u8,
}

pub fn compute_theta(
    ckey: &CommitmentKey,
    commit: &G1Projective,
    beta: &Vec<ScalarField>,
    y: ScalarField,
) -> PairingOutput<Bls12_381> {
    let n = beta.len();
    let linear_comb: G2Projective = beta
        .iter()
        .zip(ckey.u2.iter().rev())
        .map(|(&b, &u)| u * b)
        .sum();

    let pairing1 = Bls12_381::pairing(commit, linear_comb);
    let pairing2 = Bls12_381::pairing(&ckey.u1[0], &(ckey.u2[n - 1])) * y;
    pairing1 - pairing2
}

pub fn compute_M() -> G2Projective {
    G2Projective::generator()
}

pub fn gen_hash_keys() -> (ScalarField, G2Projective) {
    let mut rng = ark_std::test_rng();
    let hash_key = ScalarField::rand(&mut rng);
    let M = compute_M();
    let proj_key = M * hash_key;
    (hash_key, proj_key)
}

pub fn encrypt(
    ckey: &CommitmentKey,
    commit: &G1Projective,
    beta: &Vec<ScalarField>,
    y: ScalarField,
    message: u8,
) -> Ciphertext {
    let theta = compute_theta(ckey, commit, beta, y);
    let (hash_key, proj_key) = gen_hash_keys();
    let hash = theta * hash_key;
    let mut hash_serial: Vec<u8> = vec![];
    hash.serialize_compressed(&mut hash_serial);
    println!("hash_serial = {:?}", hash_serial);
    let l = hash.serialized_size(Compress::Yes);

    let mut rng = thread_rng();
    let rand_bytes = (0..l).map(|_| rng.gen::<u8>()).collect::<Vec<u8>>();
    let ciphertext = message
        ^ hash_serial
            .iter()
            .zip(rand_bytes.iter())
            .map(|(&a, &b)| a ^ b)
            .fold(0, |acc, x| acc ^ x);

    Ciphertext {
        proj_key,
        rand_bytes,
        ciphertext,
    }
}

pub fn decrypt(_ckey: &CommitmentKey, ct: &Ciphertext, opening: &G1Projective) -> u8 {
    let proj_hash = Bls12_381::pairing(opening, &ct.proj_key);
    let mut proj_hash_serial: Vec<u8> = vec![];
    proj_hash.serialize_compressed(&mut proj_hash_serial);
    println!("proj_hash_serial = {:?}", proj_hash_serial);
    let decrypted = ct.ciphertext
        ^ proj_hash_serial
            .iter()
            .zip(ct.rand_bytes.iter())
            .map(|(&a, &b)| a ^ b)
            .fold(0, |acc, x| acc ^ x);

    decrypted
}

#[cfg(test)]
mod test {
    use super::*;

    fn encryption_decryption(should_succeed: bool, message: u8) -> bool {
        let witness_length = 2;
        let ckey = setup_unsafe(witness_length);
        let x = vec![ScalarField::from(1), ScalarField::from(2)];
        let beta = vec![ScalarField::from(3), ScalarField::from(2)];
        let mut y = compute_func(&x, &beta);
        let (commit, r_commit) = commit(&ckey, &x);
        // println!("x = {:?}", x);
        // println!("beta = {:?}", beta);
        // println!("y = {:?}", y);
        let ct = encrypt(&ckey, &commit, &beta, y, message);
        if !should_succeed {
            y = y * y;
        }

        let opening = open(&ckey, &x, r_commit, &beta);
        let decrypted = decrypt(&ckey, &ct, &opening);
        decrypted == message
    }

    #[test]
    fn encryption_decryption_success() {
        let message = 3;
        for _ in 0..10 {
            // we're only encrypting a bit, which can be easily randomly correct, so to properly test we should run this often
            assert!(encryption_decryption(true, message));
        }
        // assert!(encryption_decryption(true, message));
    }
}
