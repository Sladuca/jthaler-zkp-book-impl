extern crate nalgebra as na;

use ark_ff::PrimeField;
use rand::prelude::*;
use na::{SMatrix, SVector};

type FieldMatrix<F, const N: usize> = SMatrix<F, N, N>;
type FieldColVec<F, const N: usize> = SVector<F, N>;

pub fn random_r<F: PrimeField>(rng: &mut ThreadRng) -> F {
    loop {
        let mut bytes = [0; 32];
        rng.fill_bytes(&mut bytes);
        let r = F::from_random_bytes(&bytes);

        if r.is_some() {
            break r.unwrap()
        }
    }
}

pub fn generate_fingerprint_vector<F: PrimeField, const N: usize>(r: F) -> FieldColVec<F, N> {
    let mut x = [F::zero(); N];

    x[0] = F::one();
    for i in 1..N {
        x[i] = x[i - 1] * r;
    }

    FieldColVec::from_row_slice(&x)
}

pub fn verify<F: PrimeField, const N: usize>(x: FieldColVec<F, N>, a: FieldMatrix<F, N>, b: FieldMatrix<F, N>, c: FieldMatrix<F, N>) -> bool {
    let y: FieldColVec<F, N> = c * x;

    let mut z = b * x;
    z = a * z;

    z == y
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;


    fn random_vec<F: PrimeField>(rng: &mut ThreadRng, len: usize) -> Vec<F> {
        let mut v = Vec::new();
        while v.len() < len {
            let mut bytes = [0; 32];
            rng.fill_bytes(&mut bytes); 
            
            match F::from_random_bytes(&bytes) {
                Some(r) => v.push(r),
                None => {}
            }
        }

        v
    }

    fn random_matrix<F: PrimeField, const N: usize>(rng: &mut ThreadRng) -> FieldMatrix<F, N> {
       let mut cols = Vec::new(); 
       for _ in 0..N {
           cols.push(FieldColVec::<F, N>::from_column_slice(&random_vec::<F>(rng, N)));
       }

       FieldMatrix::from_columns(cols.as_slice())
    }
    
    #[test]
    fn it_works() {
        const N: usize = 32;
        let mut rng = thread_rng();
        let r = random_r::<Fr>(&mut rng);
        let fingerprint_vector = generate_fingerprint_vector(r);

        let a = random_matrix::<Fr, N>(&mut rng);
        let b = random_matrix::<Fr, N>(&mut rng);

        let c = a * b;
        
        assert!(verify(fingerprint_vector, a, b, c));

        for i in 0..100 {
            let c = random_matrix::<Fr, N>(&mut rng);
            assert!(!verify(fingerprint_vector, a, b, c));
        }
    }
}
