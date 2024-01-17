use ark_bls12_381::Fr as ScalarField;
use ark_bls12_381::{Bls12_381, Fq12, G1Projective, G2Projective};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::Group;
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError};
use js_sys::Uint8Array;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::encrypt;
use crate::linear_fc;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

fn copy_vec_to_u8arr(v: &Vec<u8>) -> Uint8Array {
    let u8_arr = Uint8Array::new_with_length(v.len() as u32);
    u8_arr.copy_from(v);
    u8_arr
}

#[derive(Serialize, Deserialize)]
struct Commitment {
    commit: Vec<u8>,
    r_commit: Vec<u8>,
}

#[wasm_bindgen]
pub fn commit(ckey_bytes_1: &[u8], ckey_bytes_2: &[u8], x: u32) -> JsValue {
    let u1 = G1Projective::deserialize_compressed(ckey_bytes_1)
        .expect("deserialization should not fail");
    let u2 = G2Projective::deserialize_compressed(ckey_bytes_2)
        .expect("deserialization should not fail");
    let ckey = linear_fc::CommitmentKey {
        u1: vec![u1],
        u2: vec![u2],
    };

    let x = ScalarField::from(x);
    let (commit, r_commit) = linear_fc::commit(&ckey, &vec![x]);
    let mut commit_serial: Vec<u8> = vec![];
    commit
        .serialize_compressed(&mut commit_serial)
        .expect("serialization should not fail");
    let mut r_commit_serial: Vec<u8> = vec![];
    r_commit
        .serialize_compressed(&mut r_commit_serial)
        .expect("serialization should not fail");

    serde_wasm_bindgen::to_value(&Commitment {
        commit: commit_serial,
        r_commit: r_commit_serial,
    })
    .unwrap()
}

#[wasm_bindgen]
pub fn encrypt() {}

#[wasm_bindgen]
pub fn decrypt() {}
