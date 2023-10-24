use ndarray::{*,linalg::kron};
use ndarray_linalg::{c64,*};
use rayon::{iter::{IntoParallelRefMutIterator,ParallelIterator}, prelude::IndexedParallelIterator};
use indicatif::{ParallelProgressIterator, ProgressIterator};

use crate::{parameters::*, rad_hamiltonian::*, pa_hamiltonian::*};

/// Dispatcher that solves the TISE for the system parallelized over many different k-points
/// for a given coupling strength and command line args
pub fn rayon_dispatch(mut data: Array2<f64>, mut data_color: Array2<f64>, args:&Vec<String>, 
    k_points: &Array1<f64>, g_wc: f64, absorb_wc: Option<f64>) -> (Array2<f64>, Array2<f64>) {
    
    let prm = get_parameters(args);

    data.slice_mut(s![..,0]).assign(&k_points);

    let mut data_vec: Vec<(Array1<f64>,Array1<f64>)> = vec![(Array1::zeros(prm.nf * prm.n_kappa + 1), 
        Array1::zeros(prm.nf * prm.n_kappa + 1)); prm.nk];
    
    data_vec.par_iter_mut()
    .progress()
    .enumerate().for_each( |(k, col)| {
    // data_vec.par_iter_mut().enumerate().for_each( |(k, col)| {
        
        let mut prm_k = get_parameters(&args);

        prm_k.k = k_points[k];
        prm_k.g_wc = g_wc;

        match absorb_wc {
            Some(wc) => prm_k.wc = wc,
            None => prm_k.wc = (prm_k.wc_norm.powi(2) + (k_points[k] - prm.k_shift).powi(2)).sqrt(),
        }

        (prm_k.omega, prm_k.xi_g) = get_couplings(&prm_k);

        *col = solve_h(&prm_k);
    
    });

    for ijk in 0 .. prm.nk {
        data.slice_mut(s![ijk,1..]).assign( &(data_vec[ijk].0) );
        data_color.slice_mut(s![ijk,..]).assign( &(data_vec[ijk].1) );
    }

    (data, data_color)
}

/// Dispatcher that solves the TISE for the system iterating serially over many different k-points
/// for a given coupling strength and command line args
pub fn basic_dispatch(mut data: Array2<f64>, mut data_color: Array2<f64>, args:&Vec<String>, k_points: &Array1<f64>, g_wc: f64) -> (Array2<f64>, Array2<f64>){

    let mut prm = get_parameters(args);

    prm.g_wc = g_wc;

    data.slice_mut(s![..,0]).assign(k_points);
                
    for k in k_points.iter().progress().enumerate(){
        prm.k = *k.1;
        prm.wc = (prm.wc_norm.powi(2) + (k.1 - prm.k_shift).powi(2)).sqrt();

        (prm.omega, prm.xi_g) = get_couplings(&prm);

        let (eig_e,n_pa) = solve_h(&prm);
        data.slice_mut(s![k.0,1..]).assign(&eig_e);
        data_color.slice_mut(s![k.0,..]).assign(&n_pa);
    }

    (data, data_color)
}

/// For a given `prm: &Parameters`, this calls the function to build the Hamiltonian and then
/// solves it's eigen relations. Returns the eigen-energies and p.A photon numbers as a tuple
pub fn solve_h (prm: &Parameters) -> (Array1<f64>, Array1<f64>){

    let h = construct_h_total(prm);

    let (eig_e, eig_v) = h.eigh(UPLO::Upper).unwrap();

    let n_pa_diag = ave_photon_pa(&prm, eig_v);

    (eig_e,n_pa_diag)
}

/// Wrapper function that calls other function to generate each component of the Hamiltonian:
/// 
/// H_RAD = H_ph + T_RAD + V_RAD
/// H_C = H_ph + T_C + V 
/// 
/// Returns the total Hamiltonian as an `Array2<c64>`
pub fn construct_h_total(prm:&Parameters) -> Array2<c64> {

    match prm.hamiltonian.as_str() {
        "RAD" => return construct_rad_h(prm),
        "pA"  => return construct_pa_h(prm),
        _     => panic!("Invalid Hamiltonian type inputted")
    };
}

/// Returns the Coulomb gauge <a^+ a> for a given eigenvector (`eig_v`) in the RAD representation
/// calculated for the parameters in `prm`
fn ave_photon_pa(prm: &Parameters, eig_v:Array2<c64>) -> Array1<f64> {

    let i_m = iden(prm.n_kappa);
    let a: Array2<c64>;

    match prm.hamiltonian.as_str() {
        "RAD" => {
            a = get_a_rad(prm.nf, prm);
        },
        "pA"  => {
            a = get_a_pa(prm.nf);
        },
        _     => panic!("Invalid Hamiltonian type inputted")
    };

    let n_pa = kron(&i_m, &(&a.t()).dot(&a));
    let n_photons:Array2<c64> = (&eig_v).t().dot(&(&n_pa).dot(&eig_v));
    let n_ph_diag:Array1<f64> = n_photons.into_diag().iter().map(|&x| x.abs()).collect();

    return n_ph_diag
}

/// DEPRECIATED -- is incorrect! (I think)
/// Dispatcher that solves the TISE for the system parallelized over many different k-points
/// for a given coupling strength and command line args
pub fn _absorb_dispatch(mut data: Array2<f64>, mut data_color: Array2<f64>, args:&Vec<String>, 
    k_points: &Array1<f64>, g_wc: f64, k_ph: &f64) -> (Array2<f64>, Array2<f64>) {
    
    let prm = get_parameters(args);

    data.slice_mut(s![..,0]).assign(&k_points);

    let mut data_vec: Vec<(Array1<f64>,Array1<f64>)> = vec![(Array1::zeros(prm.nf * prm.n_kappa + 1), 
        Array1::zeros(prm.nf * prm.n_kappa + 1)); prm.nk];
    
    data_vec.par_iter_mut()
    .enumerate().for_each( |(k, col)| {
        
        let mut prm_k = get_parameters(&args);

        prm_k.g_wc = g_wc;
        prm_k.k_shift = k_points[k];
        prm_k.k = *k_ph + prm.k_shift;

        prm_k.k_ph = *k_ph;
        prm_k.wc = (prm_k.wc_norm.powi(2) + (k_ph).powi(2)).sqrt();
        
        (prm_k.omega, prm_k.xi_g) = get_couplings(&prm_k);

        *col = solve_h(&prm_k);

        let e_min = col.0[0];

        col.0 = &col.0 - e_min;
    
    });

    for ijk in 0 .. prm.nk {
        data.slice_mut(s![ijk,1..]).assign( &(data_vec[ijk].0) );
        data_color.slice_mut(s![ijk,..]).assign( &(data_vec[ijk].1) );
    }

    (data, data_color)
}