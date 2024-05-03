use blstrs::{Bls12, G1Projective, G2Projective, Gt, Scalar};
use group::{ff::Field, Curve, Group};
use pairing::Engine;

pub const L: u32 = 256; // Length of hashkey and random bits used to encrypt the message

// TODO: CommitterKey.0[n+1] is the TrapdoorKey, which should be separate.
pub type CommitterKey = (Vec<G1Projective>, Vec<G2Projective>);
pub type FunctionVars = Vec<Scalar>;
pub type FunctionCoeffs = Vec<Scalar>;
pub type FunctionOutput = Scalar;
pub type Commitment = G1Projective;

/// Generate a randomly instantiated linear functional commitment scheme.
/// Returns:
/// - CommitterKey
/// - private FunctionVars and public FunctionCoeffs, which represent the linear function
/// - FunctionOutput, which is the result of \sum_i{FunctionVars[i]*FunctionCoeffs[i]}
/// - the random scalar used to generate the CommitterKey
pub fn setup_random(
    witness_length: u64,
) -> (
    CommitterKey,
    FunctionVars,
    FunctionCoeffs,
    FunctionOutput,
    Scalar,
) {
    let mut rng = rand::thread_rng();

    // Compute scalars for the CommitterKey
    let u = Scalar::random(&mut rng);
    let mut scalars = Vec::with_capacity(2 * witness_length as usize);
    for i in 1..(2 * witness_length) + 1 {
        scalars.push(u.pow_vartime(&[i]))
    }

    // Create CommitterKey
    // TODO: u1[n+1] is the TrapdoorKey, which should be separate.
    let u1 = (0..2 * witness_length as usize)
        .into_iter()
        .map(|i| G1Projective::generator() * scalars[i])
        .collect();
    let u2 = (0..witness_length as usize)
        .into_iter()
        .map(|i| G2Projective::generator() * scalars[i])
        .collect();

    // TODO: should just do random sampling.
    let mut alpha = vec![]; // [1, 2, 1, 2, ...]
    let mut beta = vec![]; // [2, 1, 2, 1, ...]
    let mut y = Scalar::zero(); //
    for i in 0..witness_length {
        if i % 2 == 0 {
            alpha.push(Scalar::one());
            beta.push(Scalar::one().double());
        } else {
            alpha.push(Scalar::one().double());
            beta.push(Scalar::one());
        }
        y += alpha[i as usize] * beta[i as usize];
    }

    ((u1, u2), alpha, beta, y, u)
}

pub fn commit(ck: CommitterKey, alpha: FunctionVars) -> (Commitment, (FunctionVars, Scalar)) {
    let mut rng = rand::thread_rng();
    let r = Scalar::random(&mut rng);
    let r1_gen = G1Projective::generator() * r;
    let mut sum: G1Projective = G1Projective::identity();
    for i in 0..alpha.len() {
        sum += ck.0[i] * alpha[i];
    }
    let cm = sum + r1_gen;
    (cm, (alpha, r))
}

pub fn open(
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    d: (FunctionVars, Scalar),
    beta: FunctionCoeffs,
) -> G1Projective {
    let n = beta.len();
    let mut sum = G1Projective::identity();
    for i in 0..n {
        let mut w_i = ck.0[n - 1 - i] * d.1;
        for j in 0..n {
            if j != i {
                w_i += ck.0[n - i + j] * d.0[j];
            }
        }
        sum += w_i * beta[i];
    }
    sum
}

pub fn dl_of_opening(u: Scalar, d: (FunctionVars, Scalar), beta: FunctionCoeffs) -> Scalar {
    let n = beta.len();
    let mut sum = Scalar::zero();
    for i in 0..n {
        let mut dl_of_w_i = u.pow_vartime(&[(n - i) as u64]) * d.1;
        for j in 0..n {
            if j != i {
                dl_of_w_i += u.pow_vartime(&[(n + 1 - i + j) as u64]) * d.0[j];
            }
        }
        sum += dl_of_w_i * beta[i];
    }
    sum
}

pub fn dl_of_m(u: Scalar, n: u64, y: Scalar) -> Scalar {
    u.pow_vartime(&[(n + 1) as u64]) * y
}

pub fn lefthandside(
    cm: Commitment,
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    beta: FunctionCoeffs,
    n: u64,
) -> Gt {
    let mut sum = G2Projective::identity();
    for i in 0..n {
        sum += ck.1[(n - 1 - i) as usize] * beta[i as usize];
    }
    Bls12::pairing(&cm.to_affine(), &sum.to_affine())
}

pub fn righthandside(ck: (Vec<G1Projective>, Vec<G2Projective>), y: Scalar, n: u64) -> Gt {
    Bls12::pairing(
        &ck.0[0].to_affine(),
        &(ck.1[(n - 1) as usize].to_affine() * y).to_affine(),
    )
}

pub fn verify(
    ck: (Vec<G1Projective>, Vec<G2Projective>),
    cm: Commitment,
    op: G1Projective,
    beta: FunctionCoeffs,
    y: Scalar,
    n: u64,
) -> bool {
    let lefthandside = lefthandside(cm, ck.clone(), beta, n);
    let righthandside = righthandside(ck.clone(), y, n);
    let opening_pairing = Bls12::pairing(&op.to_affine(), &G2Projective::generator().to_affine());

    lefthandside == opening_pairing + righthandside
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn discreet_logs() {
        for _ in 0..128 {
            let n = 2;
            let (ck, alpha, beta, y, u) = setup_random(n);
            let (_cm, d) = commit(ck.clone(), alpha.clone());
            let op = open(ck.clone(), d.clone(), beta.clone());
            let dl_of_opening = dl_of_opening(u, d, beta.clone());
            assert_eq!(G1Projective::generator() * dl_of_opening, op);

            let righthandside = righthandside(ck.clone(), y, n);
            let dl_of_m = dl_of_m(u, n, y);
            assert_eq!(Gt::generator() * dl_of_m, righthandside);

            let opening_pairing =
                Bls12::pairing(&op.to_affine(), &G2Projective::generator().to_affine());
            let op_t = Gt::generator() * dl_of_opening;
            assert_eq!(opening_pairing, op_t);

            let sum_pairing = opening_pairing + righthandside;
            let sum_dl = dl_of_opening + dl_of_m;
            assert_eq!(Gt::generator() * sum_dl, sum_pairing);
        }
    }

    fn do_test_fc(should_succeed: bool) -> bool {
        let n = 2;
        let (ck, alpha, beta, y, _u) = setup_random(n);
        let (cm, d) = commit(ck.clone(), alpha.clone());
        let mut op = open(ck.clone(), d, beta.clone());

        if !should_succeed {
            let rng = rand::thread_rng();
            op = op * Scalar::random(rng);
        }

        verify(ck.clone(), cm, op, beta.clone(), y, n)
    }

    #[test]
    fn fc_success() {
        for _ in 0..128 {
            assert_eq!(do_test_fc(true), true);
        }
    }

    #[test]
    fn fc_failure() {
        for _ in 0..128 {
            assert_eq!(do_test_fc(false), false);
        }
    }
}
