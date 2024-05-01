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

pub fn gen_Γ_conjunction(Γ: Vec<DMatrix<BlsElement>>) -> DMatrix<BlsElement> {
    if Γ.len() != 2 {
        panic!("Currently only conjunctions of two Γ is supported");
    }

    let gt0 = BlsElement::Gt(Gt::identity());
    let g10 = BlsElement::G1Affine(G1Projective::identity().to_affine());
    let g20 = BlsElement::G2Affine(G2Projective::identity().to_affine());
    let rows = 4 * Γ.len();
    let columns = 3 * Γ.len();

    DMatrix::from_row_slice(
        rows,
        columns,
        &[
            Γ[0][0], Γ[0][4], Γ[0][8], gt0, g10, g20, Γ[0][1], Γ[0][5], Γ[0][9], gt0, g10, g20,
            Γ[0][2], Γ[0][6], Γ[0][10], gt0, g10, g20, Γ[0][3], Γ[0][7], Γ[0][11], gt0, g10, g20,
            gt0, g10, g20, Γ[1][0], Γ[1][4], Γ[1][8], gt0, g10, g20, Γ[1][1], Γ[1][5], Γ[1][9],
            gt0, g10, g20, Γ[1][2], Γ[1][6], Γ[1][10], gt0, g10, g20, Γ[1][3], Γ[1][7], Γ[1][11],
        ],
    )
}

pub fn gen_Γ_disjunction(
    Γ: Vec<DMatrix<BlsElement>>,
    θ: Vec<DMatrix<BlsElement>>,
) -> DMatrix<BlsElement> {
    if Γ.len() != 2 {
        panic!("Currently only disjunctions of two Γ is supported");
    }

    let gt0 = BlsElement::Gt(Gt::identity());
    let g10 = BlsElement::G1Affine(G1Projective::identity().to_affine());
    let g20 = BlsElement::G2Affine(G2Projective::identity().to_affine());
    let rows = 1 + 4 * Γ.len();
    let columns = 2 + 3 * Γ.len();

    DMatrix::from_row_slice(
        rows,
        columns,
        &[
            gt0, g10, g20, gt0, gt0, g10, g20, gt0, Γ[0][0], Γ[0][4], Γ[0][8], θ[0][0], gt0, g10,
            g20, gt0, Γ[0][1], Γ[0][5], Γ[0][9], θ[0][1], gt0, g10, g20, gt0, Γ[0][2], Γ[0][6],
            Γ[0][10], θ[0][2], gt0, g10, g20, gt0, Γ[0][3], Γ[0][7], Γ[0][11], θ[0][3], gt0, g10,
            g20, gt0, gt0, g10, g20, gt0, Γ[1][0], Γ[1][4], Γ[1][8], θ[1][0], gt0, g10, g20, gt0,
            Γ[1][1], Γ[1][5], Γ[1][9], θ[1][1], gt0, g10, g20, gt0, Γ[1][2], Γ[1][6], Γ[1][10],
            θ[1][2], gt0, g10, g20, gt0, Γ[1][3], Γ[1][7], Γ[1][11], θ[1][3],
        ],
    )
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

pub fn define_language_random() -> (
    (G1Affine, G2Affine, G1Affine, G2Affine),
    DMatrix<BlsElement>,
) {
    let mut rng = rand::thread_rng();
    let g1 = G1Projective::random(&mut rng).to_affine();
    let g2 = G2Projective::random(&mut rng).to_affine();
    let h1 = G1Projective::random(&mut rng).to_affine();
    let h2 = G2Projective::random(&mut rng).to_affine();

    let g1_1 = BlsElement::G1Affine(G1Affine::identity());
    let g2_1 = BlsElement::G2Affine(G2Affine::identity());
    let gt_1 = BlsElement::Gt(Gt::identity());
    let g1a = BlsElement::G1Affine(g1);
    let g2a = BlsElement::G2Affine(g2);
    let h1a = BlsElement::G1Affine(h1);
    let h2a = BlsElement::G2Affine(h2);
    let g1g2 = BlsElement::Gt(Bls12::pairing(&g1, &g2));
    let h1h2 = BlsElement::Gt(Bls12::pairing(&h1, &h2));

    let Γ = DMatrix::from_row_slice(
        4,
        3,
        vec![
            g1g2, g1_1, g2_1, gt_1, g1a, g2_1, gt_1, g1_1, g2a, h1h2, h1a, h2a,
        ]
        .as_slice(),
    );

    ((g1, g2, h1, h2), Γ)
    // Γ contains g1,g2,h1,h2, which is not very DRY but it is readable
}

pub fn gen_hash_keys(Γ: DMatrix<BlsElement>) -> (DMatrix<BlsElement>, DMatrix<BlsElement>) {
    let mut rng = rand::thread_rng();
    let a_n = DMatrix::from_fn(Γ.nrows(), 1, |_r, _c| {
        BlsElement::Scalar(Scalar::random(&mut rng))
    });
    let hp = mul_Γ(a_n.clone(), Γ);

    (a_n, hp)
}

fn msg_to_g1(msg: String) -> G1Affine {
    let empty: [u8; 0] = [];
    G1Projective::hash_to_curve(msg.as_bytes(), &empty, &empty).to_affine()
}

fn msg_to_g2(msg: String) -> G2Affine {
    let empty: [u8; 0] = [];
    G2Projective::hash_to_curve(msg.as_bytes(), &empty, &empty).to_affine()
}

pub fn generate_word(
    base_gen: (G1Affine, G2Affine, G1Affine, G2Affine),
    msg1: String,
    msg2: String,
) -> (
    (Scalar, Scalar, G1Affine, G2Affine),
    (G1Projective, G1Projective, G2Projective, G2Projective),
) {
    let msg1_g1 = msg_to_g1(msg1);
    let msg2_g2 = msg_to_g2(msg2);

    let mut rng = rand::thread_rng();
    let r1 = Scalar::random(&mut rng);
    let r2 = Scalar::random(&mut rng);

    let witness = (r1, r2, msg1_g1, msg2_g2);
    let word = (
        base_gen.0 * r1,
        (base_gen.2 * r1) + msg1_g1,
        base_gen.1 * r2,
        (base_gen.3 * r2) + msg2_g2,
    );

    (witness, word)
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

pub fn gen_θ(
    word: (G1Projective, G1Projective, G2Projective, G2Projective),
    msg1: String,
    msg2: String,
) -> DMatrix<BlsElement> {
    let msg1_g1 = msg_to_g1(msg1);
    let msg2_g2 = msg_to_g2(msg2);

    DMatrix::from_vec(
        4,
        1,
        vec![
            BlsElement::Gt(Bls12::pairing(&-word.0.to_affine(), &word.2.to_affine())),
            BlsElement::Gt(Bls12::pairing(&word.0.to_affine(), &word.3.to_affine())),
            BlsElement::Gt(Bls12::pairing(&word.1.to_affine(), &word.2.to_affine())),
            BlsElement::Gt(
                Bls12::pairing(&word.1.to_affine(), &word.3.to_affine())
                    - Bls12::pairing(&msg1_g1, &msg2_g2),
            ),
        ],
    )
}

pub fn gen_θ_disjunction(len: usize) -> DMatrix<BlsElement> {
    let mut θ: Vec<BlsElement> = vec![BlsElement::Gt(-Gt::identity())];
    for _ in 1..len {
        θ.push(BlsElement::Gt(Gt::identity()));
    }
    nalgebra::DMatrix::from_row_slice(len as usize, 1, &θ.as_slice())
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

pub fn gen_prover_hash(
    hp: DMatrix<BlsElement>,
    witnesses: Vec<(Scalar, Scalar, G1Affine, G2Affine)>,
    words: Vec<(G1Projective, G1Projective, G2Projective, G2Projective)>,
    witness_to_prove: Option<u32>,
) -> Gt {
    let mut λs: Vec<DMatrix<BlsElement>> = vec![];
    for i in 0..witnesses.len() {
        λs.push(DMatrix::from_row_slice(
            3,
            1,
            vec![
                BlsElement::Scalar(-witnesses[i].0) * BlsElement::Scalar(witnesses[i].1),
                BlsElement::Scalar(witnesses[i].0) * BlsElement::G2Affine(words[i].3.to_affine()),
                BlsElement::Scalar(witnesses[i].1) * BlsElement::G1Affine(words[i].1.to_affine()),
            ]
            .as_slice(),
        ));
    }
    let λ_final;
    match witness_to_prove {
        None => {
            let mut λ_final_vec: Vec<BlsElement> = vec![];
            for i in 0..witnesses.len() {
                for element in &λs[i] {
                    λ_final_vec.push(*element);
                }
            }
            λ_final = nalgebra::DMatrix::from_row_slice(
                λs.len() * λs[0].len(),
                1,
                &λ_final_vec.as_slice(),
            );
        }
        Some(i) => {
            let mut λ_final_vec: Vec<BlsElement> = vec![];
            if i == 0 {
                for element in &λs[i as usize] {
                    λ_final_vec.push(*element);
                }
                λ_final_vec.push(BlsElement::Scalar(-Scalar::one()));
                λ_final_vec.push(BlsElement::Scalar(Scalar::zero()));
                λ_final_vec.push(BlsElement::G2Affine(G2Affine::identity()));
                λ_final_vec.push(BlsElement::G1Affine(G1Affine::identity()));
                λ_final_vec.push(BlsElement::Scalar(Scalar::zero()));
            } else if i == 1 {
                λ_final_vec.push(BlsElement::Scalar(Scalar::zero()));
                λ_final_vec.push(BlsElement::G2Affine(G2Affine::identity()));
                λ_final_vec.push(BlsElement::G1Affine(G1Affine::identity()));
                λ_final_vec.push(BlsElement::Scalar(Scalar::zero()));
                for element in &λs[i as usize] {
                    λ_final_vec.push(*element);
                }
                λ_final_vec.push(BlsElement::Scalar(-Scalar::one()));
            } else {
                panic!("We only support choosing witness 0 or 1");
            }
            λ_final = nalgebra::DMatrix::from_row_slice(
                2 + 2 * λs[i as usize].len(),
                1,
                &λ_final_vec.as_slice(),
            );
        }
    }

    let res = mul_λ(hp, λ_final);
    if let BlsElement::Gt(result) = res {
        return result;
    }
    panic!("We did not get back a Gt element");
}

#[cfg(test)]
mod test {

    use super::*;

    fn do_sphf_pairing_test(should_succeed: bool, msg1: String, msg2: String) -> bool {
        let (base_gen, Γ) = define_language_random();

        let (hk, hp) = gen_hash_keys(Γ);

        let (w, mut word) = generate_word(base_gen, msg1.clone(), msg2.clone());

        if !should_succeed {
            // Derp the word if we need to fail
            word.1 = word.1 * w.0;
        }

        let θ = gen_θ(word, msg1.clone(), msg2.clone());
        let θ = nalgebra::DMatrix::from_row_slice(θ.len(), 1, &θ.as_slice());

        let ha = gen_verifier_hash(hk, θ);
        let hb = gen_prover_hash(hp, vec![w], vec![word], None);

        ha == hb
    }

    #[test]
    fn sphf_pairing_success() {
        assert_eq!(
            do_sphf_pairing_test(true, "Hello".to_string(), "World".to_string()),
            true
        );
    }

    #[test]
    fn sphf_pairing_failure() {
        assert_eq!(
            do_sphf_pairing_test(false, "Hello".to_string(), "World".to_string()),
            false
        );
    }

    fn do_sphf_conjunction_test(should_succeed: bool, msg1: String, msg2: String) -> bool {
        let (base1_gen, Γ1) = define_language_random();
        let (base2_gen, Γ2) = define_language_random();

        let Γ = gen_Γ_conjunction(vec![Γ1, Γ2]);

        let (hk, hp) = gen_hash_keys(Γ);

        let (w1, mut word1) = generate_word(base1_gen, msg1.clone(), msg2.clone());
        let (w2, word2) = generate_word(base2_gen, msg1.clone(), msg2.clone());

        if !should_succeed {
            // Derp the word if we need to fail
            word1.1 = word1.1 * w1.0;
        }

        let θ1 = gen_θ(word1, msg1.clone(), msg2.clone());
        let θ2 = gen_θ(word2, msg1.clone(), msg2.clone());
        let mut θ_conjunction_vec: Vec<BlsElement> = vec![];
        for element in &θ1 {
            θ_conjunction_vec.push(*element);
        }
        for element in &θ2 {
            θ_conjunction_vec.push(*element);
        }
        let θ_conjunction = nalgebra::DMatrix::from_row_slice(
            θ1.len() + θ2.len(),
            1,
            &θ_conjunction_vec.as_slice(),
        );

        let ha = gen_verifier_hash(hk, θ_conjunction);
        let hb = gen_prover_hash(hp, vec![w1, w2], vec![word1, word2], None);

        ha == hb
    }

    #[test]
    fn sphf_conjunction_success() {
        assert_eq!(
            do_sphf_conjunction_test(true, "Hello".to_string(), "World".to_string()),
            true
        );
    }

    #[test]
    fn sphf_conjunction_failure() {
        assert_eq!(
            do_sphf_conjunction_test(false, "Hello".to_string(), "World".to_string()),
            false
        );
    }

    fn do_sphf_disjunction_test(
        witness_to_prove: Option<u32>,
        word_to_mangle: Option<u32>,
        msg1: String,
        msg2: String,
    ) -> bool {
        let (base1_gen, Γ1) = define_language_random();
        let (base2_gen, Γ2) = define_language_random();

        let (w1, mut word1) = generate_word(base1_gen, msg1.clone(), msg2.clone());
        let (w2, mut word2) = generate_word(base2_gen, msg1.clone(), msg2.clone());

        match word_to_mangle {
            None => {}
            Some(word_to_mangle) => {
                if word_to_mangle == 0 {
                    // Derp the word if we need to fail
                    word1.1 = word1.1 * w1.0;
                } else if word_to_mangle == 1 {
                    word2.1 = word2.1 * w2.0;
                }
            }
        }

        let θ1 = gen_θ(word1, msg1.clone(), msg2.clone());
        let θ2 = gen_θ(word2, msg1.clone(), msg2.clone());

        let Γ = gen_Γ_disjunction(vec![Γ1, Γ2], vec![θ1.clone(), θ2.clone()]);

        let rows = Γ.nrows();
        let (hk, hp) = gen_hash_keys(Γ);

        let θ = gen_θ_disjunction(rows);

        let ha = gen_verifier_hash(hk, θ);
        let hb = gen_prover_hash(hp, vec![w1, w2], vec![word1, word2], witness_to_prove);

        ha == hb
    }

    #[test]
    fn sphf_disjunction_success_1() {
        assert_eq!(
            do_sphf_disjunction_test(Some(0), None, "Hello".to_string(), "World".to_string()),
            true
        );
    }

    #[test]
    fn sphf_disjunction_success_2() {
        assert_eq!(
            do_sphf_disjunction_test(Some(0), Some(1), "Hello".to_string(), "World".to_string()),
            true
        );
    }

    #[test]
    fn sphf_disjunction_success_3() {
        assert_eq!(
            do_sphf_disjunction_test(Some(1), None, "Hello".to_string(), "World".to_string()),
            true
        );
    }

    #[test]
    fn sphf_disjunction_success_4() {
        assert_eq!(
            do_sphf_disjunction_test(Some(1), Some(0), "Hello".to_string(), "World".to_string()),
            true
        );
    }

    #[test]
    fn sphf_disjunction_failure_1() {
        assert_eq!(
            do_sphf_disjunction_test(Some(0), Some(0), "Hello".to_string(), "World".to_string()),
            false
        );
    }

    #[test]
    fn sphf_disjunction_failure_2() {
        assert_eq!(
            do_sphf_disjunction_test(Some(1), Some(1), "Hello".to_string(), "World".to_string()),
            false
        );
    }
}
