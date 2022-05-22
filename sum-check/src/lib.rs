use ark_poly::{Polynomial, multivariate::{SparsePolynomial, SparseTerm, Term}, univariate::{SparsePolynomial as UniSparsePolynomial}, MVPolynomial};
use ark_ff::PrimeField;
use ark_std::{cfg_into_iter, cfg_iter};
use crossbeam::channel::{Sender, Receiver};
use rand::{prelude::ThreadRng, Rng};

pub type MultiPolynomial<F> = SparsePolynomial<F, SparseTerm>;
pub type UniPolynomial<F> = UniSparsePolynomial<F>;

fn bit_decompose<F: PrimeField>(mut n: u64, n_bits: usize) -> Vec<F> {
    let mut vals = Vec::with_capacity(n_bits);
    for _ in 0..n_bits {
        let val = n & 1;
        vals.push(F::from(val));
        n >>= 1;
    }
    vals
}

pub fn sum_over_boolean_hypercube<F: PrimeField>(g: &MultiPolynomial<F>) -> F {
    cfg_into_iter!(0u64..(1 << g.num_vars()))
        .map(|n| {
            let vals = bit_decompose(n, g.num_vars());
            g.evaluate(&vals)
        })
        .sum()
}


pub struct Prover<F: PrimeField> {
    g: MultiPolynomial<F>,
    challenges: Vec<F>
}

impl<F: PrimeField> Prover<F> {
    pub fn start(mut self, sum: F, sum_tx: Sender<F>, tx: Sender<UniPolynomial<F>>, rx: Receiver<F>) {
        sum_tx.send(sum).unwrap();
        for i in 0..self.g.num_vars() {
            println!("prover: {}", i);
            let g_i = self.partial_sum(i);
            if tx.send(g_i).is_err() {
                return;
            }

            match rx.recv() {
                Ok(challenge) => self.challenges.push(challenge),
                Err(_) => return
            }
        }
    }

    fn partial_sum(&self, i: usize) -> UniPolynomial<F> {
        let n_bits = self.g.num_vars() - i - 1;
        cfg_into_iter!(0u64..(1 << n_bits))
            .map(|n| bit_decompose(n, n_bits as usize))
            .fold(UniPolynomial::from_coefficients_vec(Vec::new()), |acc, point| acc + self.partial_eval(i, &point))
    }

    fn partial_eval(&self, i: usize, point: &Vec<F>) -> UniPolynomial<F> {
        cfg_iter!(self.g.terms)
            .map(|(coeff, term)| {
                let mut partial_eval_term: Option<SparseTerm> = None;
                let partial_eval_coeff = cfg_iter!(term).fold(F::one(), |coeff, (var, power)| {
                    let power = *power as u64;
                    match *var {
                        var if var == i => {
                            partial_eval_term = Some(SparseTerm::new(vec![(var, power as usize)]));
                            coeff
                        }
                        var if var < i => self.challenges[var].pow(&[power]) * coeff,
                        var => point[var - i].pow(&[power]) * coeff
                    }
                });
                (*coeff * partial_eval_coeff, partial_eval_term)
            })
            .fold(UniPolynomial::from_coefficients_vec(Vec::new()), |acc, (coeff, term)| {
                let curr_partial_eval = match term {
                    Some(term) => UniPolynomial::from_coefficients_vec(vec![(term.degree(), coeff)]),
                    None => UniPolynomial::from_coefficients_vec(vec![(0, coeff)])
                };
                acc + curr_partial_eval
            })
    }
}

pub fn build_var_degree_table<F: PrimeField>(g: &MultiPolynomial<F>) -> Vec<usize> {
    let mut table = Vec::with_capacity(g.num_vars());
    for var in 0..g.num_vars() {
        table.push(
            cfg_iter!(g.terms)
                .map(|(_, term)| {
                    cfg_iter!(term)
                        .filter(|(curr_var, _)| *curr_var == var)
                        .max_by_key(|(_, power)| *power)
                        .map_or(0, |(_, power)| *power)
                })
                .max()
                .unwrap_or(0)
        );
    }
    table
}

pub struct Verifier<F: PrimeField> {
    g: MultiPolynomial<F>,
    challenges: Vec<F>,
    target: F,
    rng: ThreadRng,
}

impl<F: PrimeField> Verifier<F> {
    pub fn start(mut self, sum_rx: Receiver<F>, tx: Sender<F>, rx: Receiver<UniPolynomial<F>>) -> bool {
        let var_degree_table = build_var_degree_table(&self.g);
        let sum = sum_rx.recv().unwrap();
        self.target = sum;
        for i in 0..self.g.num_vars() {
            let g_i = rx.recv().unwrap();

            if g_i.degree() > var_degree_table[i] {
                println!("0");
                return false;
            }
            if self.target != g_i.evaluate(&F::zero()) + g_i.evaluate(&F::one()) {
                println!("1");
                return false;
            }

            let challenge = random_field_element(&mut self.rng);
            self.challenges.push(challenge);
            
            self.target = g_i.evaluate(&challenge);

            if i == self.g.num_vars() - 1 && self.target != self.g.evaluate(&self.challenges) {
                println!("2");
                return false;
            } else if i < self.g.num_vars() - 1 {
                tx.send(challenge).unwrap();
            }
        }

        true
    }
}

pub fn random_field_element<F: PrimeField, R: Rng>(rng: &mut R) -> F {
    let mut r = None;
    while r.is_none() {
        let mut bytes = [0; 32];
        rng.fill_bytes(&mut bytes);
        r = F::from_random_bytes(&bytes);
    }

    r.unwrap()
}


#[cfg(test)]
mod tests {
    use ark_ff::Zero;
    use rand::{Rng, rngs::SmallRng, thread_rng, SeedableRng};
    use ark_bls12_381::Fr;
    use super::*;
    use crossbeam::{scope, channel};

    fn random_multivariate_polynomial<F: PrimeField>(rng: &mut ThreadRng, n_vars: usize, n_terms: usize, ) -> MultiPolynomial<F> {
        let rng = SmallRng::from_rng(rng).unwrap();
        let terms = (0..n_terms)
            .map(|_| {
                let vars = (0..n_vars)
                    .filter(|_| rng.clone().gen_bool(0.33))
                    .map(|i| (i, rng.clone().gen_range(1..=3)))
                    .collect::<Vec<_>>();
                SparseTerm::new(vars)
            })
            .map(|term| (random_field_element(&mut rng.clone()), term))
            .collect();

        MultiPolynomial::from_coefficients_vec(n_vars, terms)
    }

    fn run_protocol<F: PrimeField>(sum: F, prover: Prover<F>, verifier: Verifier<F>) -> bool {
        let (sum_tx, sum_rx) = channel::bounded(0);
        let (g_tx, g_rx) = channel::bounded(0);
        let (challenge_tx, challenge_rx) = channel::bounded(0);
        
        scope(|s| {
            s.spawn(|_| {
                prover.start(sum, sum_tx, g_tx, challenge_rx);
            });

            verifier.start(sum_rx, challenge_tx, g_rx)
        }).unwrap()
    }

    #[test]
    fn it_works() {
        let mut rng = thread_rng();
        let g = random_multivariate_polynomial::<Fr>(&mut rng, 4, 10);
        let prover = Prover { g: g.clone(), challenges: Vec::new() };
        let verifier = Verifier { g, rng, challenges: Vec::new(), target: Fr::zero() };

        let sum = sum_over_boolean_hypercube(&prover.g);
        let res = run_protocol(sum, prover, verifier);
        assert!(res);
    }

    #[test]
    fn it_works_with_bad_sum() {
        let mut rng = thread_rng();

        // use a random number instead of the actual sum
        let sum = random_field_element(&mut rng);

        let g = random_multivariate_polynomial::<Fr>(&mut rng, 4, 10);
        let prover = Prover { g: g.clone(), challenges: Vec::new() };
        let verifier = Verifier { g, rng, challenges: Vec::new(), target: Fr::zero() };

        let res = run_protocol(sum, prover, verifier);
        assert!(!res);
    }

}
