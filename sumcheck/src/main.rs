use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_bls12_381::Fr as Fq;
use ark_ff::Field;

// Common interfaces for multivariate and univariate polynomial
use ark_poly::polynomial::{DenseMVPolynomial, DenseUVPolynomial, Polynomial};

mod sumcheck;
use sumcheck::{MVPoly,UVPoly,Prover};
use sumcheck::{print_mvpoly, verify}; 

use std::io;


// Parse a polynomial from input 
fn parse_polynomial() -> MVPoly {

    println!("Input number of variable in your polynomial: ");
    
    let mut vlength = String::new();
    io::stdin().read_line(&mut vlength).expect("Failed to read line");
    let num_var: usize = vlength.trim().parse().expect("Invalid input, please enter a valid u32");

    println!("\nInput your polynomial with X_i as variables:");

    let mut input = String::new();
    let mut coef_term: Vec<(Fq, SparseTerm)> = Vec::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");
    input = input.trim().to_string();
        
    // Split terms by +
    let terms: Vec<&str> = input.split(" + ").collect();

    // Parse each term
    for term in terms{
        let parts: Vec<&str> = term.split(" * ").collect(); 

        let mut coef: u32 = 1;
        let mut spare_term_vec: Vec<(usize,usize)> = Vec::new();

        // Parse coefficient
        if let Ok(parsed_c) = parts[0].parse::<u32>(){
            coef = parsed_c;
            
            // Constant term 
            if parts.len() < 2{
                let spterm = SparseTerm::new(vec![(0, 0)]);
                coef_term.push((Fq::from(coef), spterm));
                continue;
            }
        }
        
        // parse non-constant term if any
        for p in parts{
            let vars = p;
        
            // check for non-constant (variable) part
            if let Err(parsed_p) = vars.parse::<u32>(){
            
                // Parse exponent 
                let var_expo: Vec<&str> = vars.split(" ^ ").collect();
                let mut expo: usize = 1;
                if var_expo.len() > 1{
                    expo = var_expo[1].parse().unwrap_or(1);
                }

                // Parse variable 
                let var_subscript: Vec<&str> = var_expo[0].split("_").collect();
                let subscript: usize = var_subscript[1].parse().unwrap();
                spare_term_vec.push((subscript - 1, expo));
            }

        }

        let spare_term = SparseTerm::new(spare_term_vec);

        coef_term.push((Fq::from(coef), spare_term));
    }
    
    let g = SparsePolynomial::from_coefficients_vec(num_var, coef_term);

    println!("The polynomial you input is: ");
    print_mvpoly(&g);

    g 
}



// Main function here 
fn main() {

    // Show an example 
    println!("================================================");
    println!("                Sumcheck Protocol               ");
    println!("================================================");
    println!(" ");

    println!("Example:");

    let g_example = SparsePolynomial::from_coefficients_vec(3, vec![
        (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
        (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
        (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
    ]);
    print_mvpoly(&g_example);

    let mut result = verify(&g_example, Fq::from(12));
    println!("\n --- Sumcheck Protocl Finished ---");
    println!("=> The result is {}", result); 
    println!(" ");

    let g = parse_polynomial();

    // Input for claimed sum
    println!("\nInput the claimed sum: ");
    let mut sum_input = String::new();
    io::stdin().read_line(&mut sum_input).expect("Failed to read line");
    let claimed_sum: u32 = sum_input.trim().parse().expect("Invalid input, please enter a valid u32");

    result = verify(&g, Fq::from(claimed_sum));
    println!("\n --- Sumcheck Protocl Finished ---");
    println!("=> The result is {}", result); 
    println!(" ");
}









