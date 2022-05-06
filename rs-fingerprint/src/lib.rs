use ark_ff::PrimeField;
use rand::prelude::*;

pub fn fingerprint<F: PrimeField>(file: Vec<F>, r: F) -> F {
    if file.len() == 0 {
        return F::zero();
    }

    let mut result = *file.last().unwrap();
    for &coeff in file.iter().rev().skip(1) {
        result = result * r + coeff;
    }

    result
}

pub fn random_r<F: PrimeField>(rng: &mut ThreadRng) -> F {
    let mut r = None;
    loop {
        let mut bytes = [0; 32];
        rng.fill_bytes(&mut bytes);
        r = F::from_random_bytes(&bytes);

        if r.is_some() {
            break;
        }
    }

    r.unwrap()
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;

    fn random_file<F: PrimeField>(rng: &mut ThreadRng, len: usize) -> Vec<F> {
        let mut file = Vec::new();
        while file.len() < len {
            let mut bytes = [0; 32];
            rng.fill_bytes(&mut bytes); 
            
            match F::from_random_bytes(&bytes) {
                Some(r) => file.push(r),
                None => {}
            }
        }

        file
    }

    #[test]
    fn it_works() {
        let mut rng = thread_rng();

        for _ in 0..100 {
            let alice_file = random_file::<Fr>(&mut rng, 1000);
            let bob_file = random_file::<Fr>(&mut rng, 1000);
            
            let alice_r = random_r(&mut rng);
            let alice_fingerprint = fingerprint(alice_file, alice_r);

            let bob_fingerprint = fingerprint(bob_file, alice_r);

            assert!(alice_fingerprint != bob_fingerprint);
        }
    }
}
