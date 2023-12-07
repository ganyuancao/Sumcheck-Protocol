use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_bls12_381::Fr as Fq;
use ark_ff::Field;

// Common interfaces for multivariate and univariate polynomial
use ark_poly::polynomial::{DenseMVPolynomial, DenseUVPolynomial, Polynomial};

mod sumcheck;
use sumcheck::{MVPoly,UVPoly,Prover};
use sumcheck::{print_mvpoly, verify}; 

fn main() {
    let g_example = SparsePolynomial::from_coefficients_vec(3, vec![
        (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
        (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
        (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
    ]);
    print_mvpoly(&g_example);

    let mut result = verify(&g_example, Fq::from(12));
    print!("{}", result); 


}









