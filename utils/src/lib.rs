use ark_ff::PrimeField;

pub fn bit_decompose<F: PrimeField>(mut n: u64, n_bits: usize) -> Vec<F> {
    let mut vals = Vec::with_capacity(n_bits);
    for _ in 0..n_bits {
        let val = n & 1;
        vals.push(F::from(val));
        n >>= 1;
    }
    vals
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
