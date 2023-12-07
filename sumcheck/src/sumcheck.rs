// Select a field F_q to work on (coefficient and random in F_q)
use ark_bls12_381::Fr as Fq;
use ark_ff::Field;

// Use Mutivariate and Univariate Polynomial package 
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;

// Common interfaces for multivariate and univariate polynomial
use ark_poly::polynomial::{DenseMVPolynomial, DenseUVPolynomial, Polynomial};

use ark_std::cfg_into_iter;
use rand::Rng;

// rename for Uni and Multi polynomial for convenience 
pub type MVPoly = SparsePolynomial<Fq, SparseTerm>;
pub type UVPoly = UniSparsePolynomial<Fq>;


// Converts i into an index in {0,1}^v
pub fn n_to_vec(i: usize, n: usize) -> Vec<Fq> {
	format!("{:0>width$}", format!("{:b}", i), width = n)
		.chars()
		.map(|x| if x == '1' { 1.into() } else { 0.into() })
		.collect()
}


// Construct Prover 
#[derive(Debug, Clone)]
pub struct Prover {
	pub g: MVPoly,
	pub r_vec: Vec<Fq>,
}


impl Prover {

    // The prover takes a polynomial and initalize a empty vector for randoms 
	pub fn new(g: &MVPoly) -> Self {
		Prover {
			g: g.clone(),
			r_vec: vec![],
		}
	}

    // generate a univarite polynomial 
    // g_j(X_j) = g(r_1, r_2, ..., X_j, x_{j+1}, ..., x_v) 
    // with x_{j+1}, ..., x_v in {0,1}
    pub fn obtain_unipoly(&mut self, r: Option<Fq>) -> UVPoly {
        if let Some(rr) = r {
            self.r_vec.push(rr);
        }
    
        let v = self.g.num_vars() - self.r_vec.len();
        let mut sum_poly = UVPoly::from_coefficients_vec(vec![(0, Fq::from(0))]);
    
        // Evaluate x_{j+1} to x_v
        for n in 0..(2u32.pow(v as u32 - 1)) {
            sum_poly = sum_poly + self.eval_gj(n_to_vec(n as usize, v));
        }
    
        sum_poly
    }

    // Evaluate sum on x_{j+1}, ..., x_v in {0,1}
    pub fn eval_gj(&self, points: Vec<Fq>) -> UVPoly {

         // Initialize an empty polynomial as the sum.
        let mut sum_poly = UVPoly::from_coefficients_vec(vec![]);
    
        // Iterate over each term in the polynomial g.
        for (coeff, term) in self.g.terms() {

            // Evaluate each term 
            let (coeff_eval, fixed_term) = self.eval_term(&term, &points);
            
            let curr = match fixed_term {
                // If there is no fixed term, create a constant term in the polynomial.
                None => UVPoly::from_coefficients_vec(vec![(0, *coeff * coeff_eval)]),

                // If there is a fixed term, create a term in the polynomial with the fixed term's degree.
                Some(ft) => UVPoly::from_coefficients_vec(vec![(ft.degree(), *coeff * coeff_eval)]),
            };
    
            // Summing up each term 
            sum_poly = sum_poly + curr;
        }
    
        sum_poly
    }

    // Evaluate each term in the multivariate polynomial 
    pub fn eval_term(&self, term: &SparseTerm, point: &Vec<Fq>,) -> (Fq, Option<SparseTerm>) {
        let mut fixed_term: Option<SparseTerm> = None;

        // Initialize the coefficient to 1.
        let mut coeff = Fq::from(1);
    
        // Iterate over each variable-exponent pair in the term.
        for (var, expo) in term.iter() {
            coeff = match *var {

                // Keep X_j as varaiable 
                j if j == self.r_vec.len() => {
                    fixed_term = Some(SparseTerm::new(vec![(j, *expo)]));
                    coeff
                }
                
                // Evaluate r_1,...,r_{j-1}
                j if j < self.r_vec.len() => self.r_vec[j].pow(&[*expo as u64]) * coeff,
                
                // Evaluate x_j,...,x_v
                _ => point[*var - self.r_vec.len()].pow(&[*expo as u64]) * coeff,
            };
        }
    
        // Return the resulting coefficient and fixed_term.
        (coeff, fixed_term)
    }

}


pub fn verify(g: &MVPoly, c_1: Fq) -> Fq{
    let mut p = Prover::new(g);
    let mut g1 = p.obtain_unipoly(None);
    let mut expected_c = g1.evaluate(&Fq::from(0)) + g1.evaluate(&Fq::from(1));
    expected_c
}



// Utility Here 
// Print multi-variate polynomial 
pub fn print_mvpoly(poly: &MVPoly) {
    // Iterate through the terms of the polynomial
    print!("Polynomial to be verified g =");
    for (i,(coeff, term)) in poly.terms().iter().enumerate() {
        print_term(coeff, term);
        if i < poly.terms().len()-1{
            print!(" +");
        }
    }
    println!(); // Print a newline at the end
}

pub fn print_term(coeff: &Fq, term: &SparseTerm) {
    // Constant Term
    if term.is_empty() {
        print!("{}", coeff);
    } 
    // Non-constant term
    else {
        for (i, (variable, exponent)) in term.iter().enumerate() {
            if *exponent == 1 {
                print!(" X_{}", variable);
            } else {
                print!(" X_{}^{}", variable, exponent);
            }

             // Print '*' only if there are more terms to come
             if i < term.len() - 1 {
                print!(" *");
            }
        }
    }
}

