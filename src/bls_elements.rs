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

// Wrapper element allows for creating matrices of a single type
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum BlsElement {
    Scalar(Scalar),
    G1Affine(G1Affine),
    G2Affine(G2Affine),
    Gt(Gt),
}

impl std::ops::Mul for BlsElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        if let BlsElement::Scalar(element1) = self {
            if let BlsElement::Scalar(element2) = rhs {
                return BlsElement::Scalar(element1 * element2);
            } else if let BlsElement::Gt(element2) = rhs {
                return BlsElement::Gt(element2 * element1);
            } else if let BlsElement::G2Affine(element2) = rhs {
                return BlsElement::G2Affine((element2 * element1).to_affine());
            } else if let BlsElement::G1Affine(element2) = rhs {
                return BlsElement::G1Affine((element2 * element1).to_affine());
            }
        }
        if let BlsElement::G1Affine(element1) = self {
            if let BlsElement::Scalar(element2) = rhs {
                return BlsElement::G1Affine((element1 * element2).to_affine());
            } else if let BlsElement::G2Affine(element2) = rhs {
                return BlsElement::Gt(Bls12::pairing(&element1, &element2));
            }
        }
        if let BlsElement::G2Affine(element1) = self {
            if let BlsElement::Scalar(element2) = rhs {
                return BlsElement::G2Affine((element1 * element2).to_affine());
            } else if let BlsElement::G1Affine(element2) = rhs {
                return BlsElement::Gt(Bls12::pairing(&element2, &element1));
            }
        }
        if let BlsElement::Gt(element1) = self {
            if let BlsElement::Scalar(element2) = rhs {
                return BlsElement::Gt(element1 * element2);
            }
        }

        panic!(
            "BlsElement types can't be multiplied - not recognized {:?} {:?}",
            self, rhs
        );
    }
}

impl std::ops::Add for BlsElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        if let BlsElement::Scalar(element1) = self {
            if let BlsElement::Scalar(element2) = rhs {
                return BlsElement::Scalar(element1 + element2);
            }
        }

        // if let BlsElement::G1Affine(element1) = self {
        //     if let BlsElement::G1Affine(element2) = rhs {
        //         let element2projective =
        //             G1Projective::from_compressed(&element2.to_compressed()).unwrap();
        //         return BlsElement::G1Affine((element1 + element2projective).to_affine());
        //     }
        // }
        // if let BlsElement::G2Affine(element1) = self {
        //     if let BlsElement::G2Affine(element2) = rhs {
        //         let element2projective =
        //             G2Projective::from_compressed(&element2.to_compressed()).unwrap();
        //         return BlsElement::G2Affine((element1 + element2projective).to_affine());
        //     }
        // }

        if let BlsElement::Gt(element1) = self {
            if let BlsElement::Gt(element2) = rhs {
                return BlsElement::Gt(element1 + element2);
            }
        }

        panic!(
            "BlsElement types can't be added - not recognized {:?} {:?}",
            self, rhs
        );
    }
}

impl std::ops::Sub for BlsElement {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        if let BlsElement::Gt(element1) = self {
            if let BlsElement::Gt(element2) = rhs {
                return BlsElement::Gt(element1 - element2);
            }
        }

        panic!(
            "BlsElement types can't be added - not recognized {:?} {:?}",
            self, rhs
        );
    }
}
