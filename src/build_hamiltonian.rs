use std::f64::consts::PI;
use ndarray::{*,linalg::kron};
use ndarray_linalg::{c64,*};
use statrs::function::gamma::gamma_ui;


use crate::parameters::*;

/// Wrapper function that calls other function to generate each component of the Hamiltonian:
/// 
/// H_RAD = H_ph + T_RAD + V_RAD
/// H_C = H_ph + T_C + V 
/// 
/// Returns the total Hamiltonian as an `Array2<c64>`
pub fn construct_h_total(prm:&Parameters) -> Array2<c64> {

    match prm.hamiltonian.as_str() {
        "RAD" => return construct_rad_h(prm),
        "pA"  => panic!("Coulomb Gauge Hamiltonian not implemented yet"),
        _     => panic!("Invalid Hamiltonian type inputted")
    };
}

fn construct_rad_h(prm:&Parameters) -> Array2<c64> {
    // Define identities: 
    let i_m = iden(prm.n_kappa);


    let k_e = get_kinetic(&prm);
    let v_shifted = get_shifted_v(&prm);
    let h_ph = get_h_ph(&prm);

    let mut h_total: Array2<c64> = Array2::zeros((prm.nf*prm.n_kappa,prm.nf*prm.n_kappa));

    h_total = h_total + kron(&i_m, &h_ph);
    h_total = h_total + k_e;
    h_total = h_total + v_shifted;

    h_total
}

/// Returns H_ph for a given set of parameters `prm`
fn get_h_ph(prm:&Parameters) -> Array2<c64>{

    let range = to_c_1(Array1::linspace(0.0, (prm.nf - 1) as f64 , prm.nf));

    let h_ph: Array2<c64> = Array2::from_diag(&range) * prm.omega;

    h_ph
}

/// Calculates T_RAD for a given set of parameters `prm`
fn get_kinetic(prm:&Parameters) -> Array2<c64>{
    
    let i_m = iden(prm.n_kappa);
    let i_ph = iden(prm.nf);

    let m_eff = c64::from(prm.m * (1.0 + 2.0 * prm.g_wc.powi(2)));

    let a = get_a(prm.nf, prm);

    let p = to_c_2(Array2::from_diag(&(prm.kappa_grid.clone() + prm.k)));

    let mut p_new = c64::from(prm.hbar) * kron(&p, &i_ph);
    p_new = p_new - kron(&i_m,  &(a.t().dot(&a))) * c64::from(prm.hbar*(prm.k - prm.k_shift));

    let k_e = &p_new.dot(&p_new) / (2.0 * m_eff);

    k_e
}

/// Calculates T_RAD for a given set of parameters `prm`
fn get_shifted_v(prm:&Parameters) -> Array2<c64>{
    
    // Generate Chi
    let b = get_b(prm.nf);
    let chi = &b.t() + &b;

    // Initialize constants as c64
    let xi_g = c64::from(prm.xi_g);
    let z = c64::from(prm.z);

    let n = prm.nf * prm.n_kappa;

    let kappa_grid = to_c_1(prm.kappa_grid.clone());
    let kappa_grid2 = to_c_1(prm.kappa_grid2.clone());

    // Generate V with phase factor
    let mut v_shifted: Array2<c64> = Array2::zeros((n,n));

    let kappa_max = kappa_grid.norm_max();

    
    for (l, k2) in kappa_grid2.iter().enumerate() {

        let phase= &chi * xi_g * kappa_grid2[l] * c64::i();

        let mut exponent = iden(prm.nf);
        if phase.map(|x| x.abs()).opnorm_one().unwrap() > f64::EPSILON * 2. {
            (exponent, _) = ndarray_linalg::expm(&phase); // This is still unverified but probably fine
        }
        
        for (m, k1) in kappa_grid.iter().enumerate(){
            let mut v_diff: Array2<c64> = Array2::zeros((prm.n_kappa,prm.n_kappa));
            if (k1 + k2).abs() <= kappa_max {
                if k2.abs() <= 1e-7 {
                    v_diff[[m,m]] = c64::from(0.0);
                    // println!("Zero!");
                }
                else {
                    let shift_ind = (prm.n_kappa2 -1) / 2;
                    v_diff[[m,m + l - shift_ind]] =  - z / 2.0 / PI * gamma_ui(1e-10, (k2.re() / 2.0 / prm.r_0).powi(2));
                    // println!("vdiff = {}",v_diff[[m,m + l - shift_ind]]);
                }
            }

            v_shifted = v_shifted + kron(&v_diff, &exponent);
        }
    }
    v_shifted
}

/// Calculates Bogoliubov transformed frequency, Ω, and effective coupling, ξ, 
/// from g_wc and wc
pub fn get_couplings(prm:&Parameters) -> (f64,f64){

    let omega = (prm.wc.powi(2) + 2.0 * prm.g_wc.powi(2)).sqrt();
    let x_omega = (prm.hbar / (prm.m * prm.omega)).sqrt();
    let xi_g = prm.g_wc * x_omega / omega;

    (omega, xi_g)
}

/// Calculates the Bogoliubov transformed b operator matrix
fn get_b(nf:usize) -> Array2<c64>{
    let mut b: Array2<c64> = Array2::zeros([nf,nf]);

    for m in 0..nf - 1 {
        let val = ((m + 1) as f64).sqrt();
        b[[m , m+1]] =c64::from(val);
    }

    b
}

/// Calculates the Coulomb gauge a operator matrix in the b Fock state basis
pub fn get_a(nf:usize, prm:&Parameters) -> Array2<c64>{
    let wc = prm.wc;
    let omega = prm.omega;

    let s_wc = (omega / wc).sqrt();

    let u = c64::from(0.5 * ( s_wc - 1.0 / s_wc ));
    let v = c64::from(0.5 * ( s_wc + 1.0 / s_wc ));

    let b = get_b(nf);

    -u * &b.t() + v * &b 
}

/// Converts an `Array1<f64>` to an `Array1<c64>`
fn to_c_1(arr:Array1<f64>) -> Array1<c64> {
    let new_arr: Array1<c64> = arr.iter().map(|&e| c64::from(e)).collect();
    new_arr
}

/// Converts an `Array2<f64>` to an `Array2<c64>`
fn to_c_2(arr:Array2<f64>) -> Array2<c64>{

    let new_arr = arr.map(|x| c64::from(x));

    new_arr    
}

/// Finds the identity matrix of size `n`
pub fn iden(n:usize) -> Array2<c64>{
    let one_d: Array1<f64> = Array1::ones(n);
    let one_d_c = to_c_1(one_d);

    let iden = Array2::from_diag(&one_d_c);
    iden
}