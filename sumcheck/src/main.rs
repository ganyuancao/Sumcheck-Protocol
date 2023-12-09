use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_bls12_381::Fr as Fq;
use ark_ff::Field;

// Common interfaces for multivariate and univariate polynomial
use ark_poly::polynomial::{DenseMVPolynomial, DenseUVPolynomial, Polynomial};

mod sumcheck;
use sumcheck::{MVPoly,UVPoly,Prover};
use sumcheck::{print_mvpoly, verify}; 

use std::io;

fn main() {

    // Show an example 
    println!("================================");
    println!("        Sumcheck Protocol       ");
    println!("================================");
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

    let g = parse_poly();

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



// Function to parse user input and extract Fq value and SparseTerm vector
fn parse_poly() -> MVPoly {

    println!("Input number of variable in your polynomial: ");
    
    let mut vlength = String::new();
    io::stdin().read_line(&mut vlength).expect("Failed to read line");
    let num_var: usize = vlength.trim().parse().expect("Invalid input, please enter a valid u32");

    
    println!("\nInput your polynomial term by term here following the format ");
    println!("coefficient, [(variable,exponent),(variable,exponent),...]");
    println!("e.g. 2, [(0,1),(1,3)] means 2 * X_1 * X_2^3");
    println!("End by inputing an empty line.");

    let mut input = String::new();
    let mut coef_term: Vec<(Fq, SparseTerm)> = Vec::new();
    
    // Read lines until an empty line is encountered
    while let Ok(_) = std::io::stdin().read_line(&mut input) {
        input = input.trim().to_string();
        
        // Break if an empty line is encountered
        if input.is_empty() {
            break;
        }

        // Parse the input and create the corresponding Fq and SparseTerm elements
        let parts: Vec<&str> = input.split(", ").collect();
        let coef: u32 = parts[0].trim().parse().expect("Error parsing coefficient.");

        let term_str = parts[1].trim_matches(|c| c == '[' || c == ']').trim();
    
        let term_vec: Vec<(usize, usize)> = term_str.split("),").map(|pair| {
            let pair_str = pair.trim_matches(|c| c == '(' || c == ')');
            let pair_parts: Vec<&str> = pair_str.split(",").collect();
            let first = pair_parts[0].trim().parse().expect("Error parsing first value of pair");
            let second = pair_parts[1].trim().parse().expect("Error parsing second value of pair");
            (first, second)
        }).collect();

        let spterm = SparseTerm::new(term_vec);

        coef_term.push((Fq::from(coef), spterm)); 
        
        input.clear();
    }

    let g = SparsePolynomial::from_coefficients_vec(num_var, coef_term);

    println!("The polynomial you input is: ");
    print_mvpoly(&g);

    g 
}









