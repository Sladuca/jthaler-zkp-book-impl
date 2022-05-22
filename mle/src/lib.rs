use ark_ff::PrimeField;
use ark_std::cfg_into_iter;

pub fn precompute_f_ws<F, Func>(f: Func, n_bits: usize) -> Vec<F>
where
    F: PrimeField,
    Func: Fn(u64) -> F
{
    if n_bits > 64 {
        panic!("n_bits must be <= 64");
    }

    let mut f_ws = Vec::new();
    let mut w: u64 = 0;
    for _ in 0..(1 << n_bits) {
        f_ws.push(f(w));
        w += 1;
    }
    
    f_ws
}

fn chi_w<F: PrimeField>(mut w: u64, x: &Vec<F>) -> F
{
    let n_bits = x.len();
    let mut prod = F::one();
    for i in 0..n_bits {
        prod *= chi_term(F::from(w & 1), x[i]);
        w >>= 1;
    }

    prod
}

fn chi_term<F: PrimeField>(w_bit: F, x: F) -> F {
    w_bit * x + (F::one() - x) * (F::one() - w_bit)
}

// assumes f_ws.len() is a power of two
pub fn eval_mle_naive<F: PrimeField>(f_ws: &[F], x: &Vec<F>) -> F
{
    f_ws.iter()
        .enumerate()
        .map(|(i, &f_w)| f_w * chi_w(i as u64, x))
        .sum()
}


pub fn build_chi_table<F: PrimeField>(f_ws: &[F], x: &Vec<F>, n_bits: usize) -> Vec<F> {
    match n_bits {
        1 => vec![chi_term(F::from(0u64), x[n_bits - 1]), chi_term(F::from(1u64), x[n_bits - 1])],
        _ => cfg_into_iter!(build_chi_table(f_ws, x, n_bits - 1))
            .flat_map(|chi_inner| [
                chi_inner * chi_term(F::from(0u64), x[n_bits - 1]),
                chi_inner * chi_term(F::from(1u64), x[n_bits - 1])
            ])
            .collect()
    }
}

pub fn eval_mle_memo<F: PrimeField>(f_ws: &[F], x: &Vec<F>) -> F {
    let chi_table = build_chi_table(f_ws, x, x.len());
    f_ws
        .iter()
        .zip(chi_table.iter())
        .map(|(&f_w, &chi_w)| f_w * chi_w)
        .sum()
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;
    use utils::bit_decompose;

    fn test_binary_fn(x: u64) -> Fr {
        Fr::from(x)
    }

    #[test]
    fn it_works() {
        const N_BITS: usize = 16;
        let fws = precompute_f_ws(&test_binary_fn, N_BITS);

        let test_vals: Vec<u64> = vec![0, 7, 69, 420, 1453];

        // test naive
        for &val in test_vals.iter() {
            let x = bit_decompose(val, N_BITS);
            let mle_eval = eval_mle_naive(&fws, &x);
            assert_eq!(mle_eval, Fr::from(val));
        }

        // test memo
        for val in test_vals {
            let mut x = bit_decompose(val, N_BITS);
            x.reverse();
            let mle_eval = eval_mle_memo(&fws, &x);
            assert_eq!(mle_eval, Fr::from(val));
        }
    }
}
