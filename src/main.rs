mod build_hamiltonian;
use build_hamiltonian::*;

pub mod parameters;
use parameters::*;

// pub mod expm;
// use expm::*;

use statrs::function::gamma::gamma_ui;

fn main() {
    ///////////// Define Parameters //////////////

    let mut prm = get_parameters();

    // let hph = get_h_ph(&prm);
    
    // let b = get_b(prm.nf);

    // (prm.omega, prm.xi_g) = get_couplings(&prm);

    // let a = get_a(prm.nf,&prm);
    // // let a: = get_a(nf, prm);
    // println!("{}", gamma_ui(1e-10, 1.0e-15));
    // println!("{}",hph);
}

