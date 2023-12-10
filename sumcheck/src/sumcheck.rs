// Select a field F_q to work on (coefficient and random in F_q)
use ark_bls12_381::Fr as Fq;
use ark_ff::Field;

// Use Mutivariate and Univariate Polynomial package 
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;

// Common interfaces for multivariate and univariate polynomial
use ark_poly::polynomial::{DenseMVPolynomial, DenseUVPolynomial, Polynomial};

use rand::Rng;

// rename for Uni and Multi polynomial for convenience 
pub type MVPoly = SparsePolynomial<Fq, SparseTerm>;
pub type UVPoly = UniSparsePolynomial<Fq>;

use std::collections::HashMap;
use std::str::Chars;


// Converts i into an vector in {0,1}^v
pub fn n_to_vec(i: usize, n: usize) -> Vec<Fq> {
    let binary_str = format!("{:0>width$}", format!("{:b}", i), width = n);

    binary_str.chars().map(|x| 
        if x == '1' { 
            Fq::from(1) 
        } 
        else { 
            Fq::from(0) 
        }).collect()
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

// Compute deg_j(g)
pub fn degj(g: &MVPoly) -> Vec<usize> {
    let mut deg: Vec<usize> = vec![0; g.num_vars()];

    for (_, term) in g.terms().iter() {
        for (var, exp) in term.iter() {
            deg[*var] = deg[*var].max(*exp);
        }
    }

    deg
}

// Get a random value in the field Fq
pub fn get_rand() -> Option<Fq> {
	let mut rng = rand::thread_rng();
	let r: Fq = rng.gen();
	Some(r)
}


// Verify steps 
pub fn verify(g: &MVPoly, c_1: Fq) -> bool{

    println!("\n --- Running Sumcheck Protocol --- \n");

    let mut prover = Prover::new(g);

    // First round
    let g_1 = prover.obtain_unipoly(None);
    let g1_sum = g_1.evaluate(&Fq::from(0)) + g_1.evaluate(&Fq::from(1));
    let deg_xj = degj(&g);

    println!("Round {}", 1); 
    println!("C_1 = {}", c_1); 
    println!("g_1(0) + g(1) = {}", g1_sum);
    
    if c_1 != g1_sum || g_1.degree() > deg_xj[0]{
        return false;
    }
    
    println!(" ");

    let v = prover.g.num_vars(); 

    let mut g_jm1 = g_1;

    // Round 2 to v-1
    for j in 1..v{
        println!("Round {}", j + 1); 
        
        let r = get_rand();
        let g_j = prover.obtain_unipoly(r);
        
        let gj_sum = g_j.evaluate(&Fq::from(0)) + g_j.evaluate(&Fq::from(1));
        let g_jm1evalr = g_jm1.evaluate(&r.unwrap());
        
        println!("r_{} = {}", j, &r.unwrap()); 
        println!("g_{}(r_{}) = {}", j, j, g_jm1evalr); 
        println!("g_{}(0) + g_{}(1) = {}", j + 1, j + 1, gj_sum); 

        if gj_sum != g_jm1evalr || g_j.degree() > deg_xj[j]{
            return false;
        }

        println!(" ");
        
        g_jm1 = g_j; 
    }

    println!("Round {}", v + 1); 
    let r = get_rand();
    // g_v(r_v) here
    let g_jm1evalr = g_jm1.evaluate(&r.unwrap());

    // g(r_1, ..., r_v) here 
    prover.r_vec.push(r.unwrap());
	let g_r1rv = prover.g.evaluate(&prover.r_vec);

    println!("r_{} = {}", v, &r.unwrap()); 
    println!("g_{}(r_{}) = {}", v, v, g_jm1evalr); 
    println!("g(r_1,...,r_{}) = {}", v, g_r1rv); 

    if g_jm1evalr != g_r1rv{
        return false;
    }

    true
}


// Print multi-variate polynomial 
pub fn print_mvpoly(poly: &MVPoly) {
    // Iterate through the terms of the polynomial
    print!("g = ");
    for (i,(coeff, term)) in poly.terms().iter().enumerate() {
        print_term(coeff, term);
        if i < poly.terms().len()-1{
            print!(" + ");
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
        if *coeff != Fq::from(1){
            print!("{} * ", coeff);
        }
        for (i, (variable, exponent)) in term.iter().enumerate() {
            if *exponent == 1 {
                print!("X_{}", variable + 1);
            } else {
                print!("X_{} ^ {}", variable + 1, exponent);
            }

             // Print '*' only if there are more terms to come
             if i < term.len() - 1 {
                print!(" * ");
            }
        }
    }
}

