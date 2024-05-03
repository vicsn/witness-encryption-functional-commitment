# Witness encryption for succinct functional commitments

This library contains a small implementation of witness encryption for succinct functional commitments for [Libert et al.’s FC](https://eprint.iacr.org/2016/766).

> :exclamation: This library has not been audited and should not be used in production

Witness encryption allows you to encrypt a value to whoever has knowledge of a witness. More concretely, in this work, the encryption statement consists of a commitment `cm`, a function `G` and a value `y`; the decryption witness consists of a (non succinct) NIZK proof about the fact that `cm` opens to `v`  such that `y=G(v)`.

The above is made possible through two main building blocks:
1. Smooth Projective Hash Functions
2. Functional Commitments

To see how it can be used, run `cargo test encryption_decryption`.

The SPHF tests use examples from Fabrice Benhamouda's thesis: Diverse modules and zero-knowledge.

## Credits
A big thank you to the authors of [Witness Encryption for Succinct Functional Commitments and Applications](https://eprint.iacr.org/2022/1510) and its sources. Moreover, this existing [SPHF repo](https://github.com/oblazy/rust-sphf/) was a great help.