// use csv::{ReaderBuilder,WriterBuilder};
// use ndarray::Array2;
// use std::{fs::File, result::Result, path::Path, env, time::Instant};
// use ndarray_csv::{Array2Writer, Array2Reader, ReadError};
use std::{env, time::Instant};

mod rad_hamiltonian;
mod pa_hamiltonian;

pub mod parameters;
use parameters::*;

mod solve_hamiltonian;
// use solve_hamiltonian::*;

mod routines;
use routines::*;

mod graphing;
// use graph_disp::*;

/*

Rust implementation of the RAD Hamiltonian beyong the long-wavelength approximation
for k modes coupled to a single particle. This program solves the dispersion 
relations of a single electron in an infinite periodic erf() potential for multiple
coupling strengths and k-points with the option of displacing the photonic k from 
the matter k using the k-shift parameter.

To run perform the following commands:

cargo build --release
./target/release/rad-in-rust {log_g_min} {log_g_max} {ng} {nk} {n_kappa} {nf}

For use on bluehive, the parallel.py script dispatches this routine accross multiple
nodes.

*/


/**
# Main function 
Reads the command line arguments to generate parameters and run the specified routine

Command line input has the form:

`routine: String`,
`wc_norm: f64
`log_g_min: f64`,
`log_g_max: f64`,
`ng: usize`,
`nk: usize`,
`n_kappa: usize`,
`nf: usize`
*/
fn main() {

    env::set_var("RUST_BACKTRACE", "1");
    
    // Start the timer to determine
    let now = Instant::now();

    // Read command line arguments: routine, log_g_min, log_g_max, ng, nk, n_kappa, nf
    let args: Vec<String> = env::args().collect();
    println!("{}", &args[1]);
    
    // Case when routine == disp
    if args[1] == "disp" {
        let args_disp = &args[1..].to_vec();

        // Initialize the prm struct based off the command-line arguments
        let prm = get_parameters(args_disp);

        get_disps(prm, args_disp);
    }
    else if args[1] == "absorb" {
        let args_abs = &args[1..].to_vec();
        
        // Initialize the prm struct based off the command-line arguments
        let prm = get_parameters(args_abs);

        get_absorption(prm, &args_abs);
    }
    else if args[1] == "disp2d" {
        let args_abs = &args[1..].to_vec();
        
        // Initialize the prm struct based off the command-line arguments
        let prm = get_parameters(args_abs);

        get_disps_2_d(prm, &args_abs);
    }
    else {
        panic!("{} is not a valid routine", args[1])
    }

    println!("Time elapsed = {} sec", now.elapsed().as_secs());
}

