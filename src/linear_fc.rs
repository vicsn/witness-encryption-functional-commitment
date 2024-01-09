use blstrs::Bls12;
use blstrs::G1Projective;
use blstrs::G2Projective;
use blstrs::Gt;
use blstrs::Scalar;

use group::ff::Field;
use group::Curve;
use group::Group;
use pairing::Engine;
use std::convert::TryInto;

pub const L: u32 = 256; // Length of hashkey and random bits used to encrypt the message

pub fn setup_random(
    witness_length: u64,
) -> (
    (Vec<G1Projective>, Vec<G2Projective>),
    Vec<Scalar>,
    Vec<Scalar>,
    Scalar,
) {
    let mut rng = rand::thread_rng();
    let u = Scalar::random(&mut rng);
    let u1_gen = G1Projective::generator() * u;
    let u2_gen = G2Projective::generator() * u;

    let mut u1 = vec![];
    let mut u2 = vec![];

    // Create a list of Scalars with values [1,2n]
    let mut scalars = vec![];
    let mut base_uint = 1;
    let mut base_scalar = Scalar::one();
    for i in 1..1 + (2 * witness_length) {
        if i == base_uint {
            scalars.push(base_scalar);
            base_uint = base_uint * 2;
            base_scalar = base_scalar.double();
        } else if i < base_uint {
            let diff: usize = (base_uint - i).try_into().unwrap();
            scalars.push(base_scalar - scalars[(diff - 1) as usize]);
        } else if i > base_uint {
            panic!("this should not be able to happen");
        }
    }

    for i in 1..1 + (2 * witness_length) {
        // NOTE: even though we don't use scalars[n+1], we add it to our scalars for developer-friendliness
        u1.push(u1_gen * scalars[(i - 1) as usize]);
    }
    for i in 1..1 + witness_length {
        u2.push(u2_gen * scalars[(i - 1) as usize]);
    }

    let mut α = vec![];
    let mut β = vec![];
    let mut y = Scalar::zero();
    for i in 0..witness_length {
        if i % 2 == 0 {
            α.push(Scalar::one());
            β.push(Scalar::one().double());
        } else {
            α.push(Scalar::one().double());
            β.push(Scalar::one());
        }
        y += α[i as usize] * β[i as usize];
    }

    ((u1, u2), α, β, y)
}

pub fn commit(
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    α: Vec<Scalar>,
) -> (G1Projective, (Vec<Scalar>, Scalar)) {
    let n = α.len();
    let mut rng = rand::thread_rng();
    let r = Scalar::random(&mut rng);
    let r1_gen = G1Projective::generator() * r;
    let mut sum: G1Projective = G1Projective::identity();
    for i in 0..n {
        sum += ck.0[i as usize] * α[i as usize];
    }
    (sum + r1_gen, (α, r))
}

pub fn open(
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    d: (Vec<Scalar>, Scalar),
    β: Vec<Scalar>,
) -> G1Projective {
    let n = β.len();
    let mut sum = G1Projective::identity();
    for i in 0..n {
        //let mut Wi = ck.0[(n-i-1) as usize]
        let mut ADJUSTMENT_1 = G1Projective::generator();
        if i == 0 {
            ADJUSTMENT_1 = ADJUSTMENT_1 * Scalar::one().double();
        }
        let mut Wi = ADJUSTMENT_1 * d.1;
        for j in 0..n {
            if j != i {
                Wi += ck.0[(n + j - i) as usize] * d.0[j as usize];
            }
        }
        sum += Wi * β[i as usize];
    }
    sum
}

pub fn verification_pairings(
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    cm: G1Projective,
    op: G1Projective,
    β: Vec<Scalar>,
    y: Scalar,
    n: u64,
    α: Vec<Scalar>,
) -> (Gt, Gt, Gt) {
    let mut sumb = G2Projective::identity();
    for i in 0..n {
        sumb += ck.1[(n - 1 - i) as usize] * β[i as usize];
    }

    let ADJUSTMENT_2 = ck.1[0].to_affine(); // G^u instead of G
    let ADJUSTMENT_3 = ck.0[0].to_affine() * β[1] * α[0]; // G^(u*α1*β2)

    let lefthandside = Bls12::pairing(&cm.to_affine(), &sumb.to_affine());
    let righthandside1 = Bls12::pairing(&(op - ADJUSTMENT_3).to_affine(), &ADJUSTMENT_2);
    let righthandside2 = Bls12::pairing(
        &ck.0[0].to_affine(),
        &(ck.1[(n - 1) as usize].to_affine() * y).to_affine(),
    );
    (lefthandside, righthandside1, righthandside2)
}

#[cfg(test)]
mod test {

    use super::*;

    fn do_test_fc(should_succeed: bool) -> bool {
        let n = 2;
        let (ck, α, β, y) = setup_random(n);
        let (cm, d) = commit(ck.clone(), α.clone());
        let mut op = open(ck.clone(), d, β.clone());

        if !should_succeed {
            let rng = rand::thread_rng();
            op = op * Scalar::random(rng);
        }

        let (e1, e2, e3) = verification_pairings(ck.clone(), cm, op, β.clone(), y, n, α.clone());
        // println!("e2: {:?}", e2);
        // println!("e2*y: {:?}", e2*y);
        let test_identity = G1Projective::identity();
        println!("test_identity.double(): {:?}", test_identity.double());
        println!(
            "test_identity*Scalar::one().double(): {:?}",
            test_identity * (Scalar::one().double())
        );

        Gt::identity() == (e1 - (e2 + e3))
        // e1 == e2 + e3
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
