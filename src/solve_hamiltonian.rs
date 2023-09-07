use ndarray::{*,linalg::kron};
use ndarray_linalg::{c64,*};
use rayon::{iter::IntoParallelRefMutIterator, prelude::IndexedParallelIterator};
use rayon::iter::ParallelIterator;
use indicatif::{ParallelProgressIterator, ProgressIterator};

use crate::{parameters::*, build_hamiltonian::*};


/// Dispatcher that solves the TISE for the system parallelized over many different k-points
/// for a given coupling strength and command line args
pub fn rayon_dispatch(mut data: Array2<f64>, mut data_color: Array2<f64>, args:&Vec<String>, 
    k_points: &Array1<f64>, g_wc: f64) -> (Array2<f64>, Array2<f64>) {
    
    let prm = get_parameters(args);

    data.slice_mut(s![..,0]).assign(&k_points);

    let mut data_vec: Vec<(Array1<f64>,Array1<f64>)> = vec![(Array1::zeros(prm.nf * prm.n_kappa + 1), 
        Array1::zeros(prm.nf * prm.n_kappa + 1)); prm.nk];
    
    data_vec.par_iter_mut().progress().enumerate().for_each( |(k, col)| {
        
        let mut prm_k = get_parameters(&args);

        prm_k.k = k_points[k];
        prm_k.g_wc = g_wc;

        prm_k.wc = (prm_k.wc_norm.powi(2) + (k_points[k] - prm_k.k_shift).powi(2)).sqrt();

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
        data_color.slice_mut(s![k.0,1..]).assign(&n_pa);
    }

    (data, data_color)
}

/// For a given `prm: &Parameters`, this calls the function to build the Hamiltonian and then
/// solves it's eigen relations. Returns the eigen-energies and p.A photon numbers as a tuple
fn solve_h (prm: &Parameters) -> (Array1<f64>, Array1<f64>){

    let h = construct_h_total(prm);

    let (eig_e, eig_v) = h.eigh(UPLO::Upper).unwrap();

    let n_pa_diag = ave_photon_pa(&prm, eig_v);

    (eig_e,n_pa_diag)
}

/// Returns the Coulomb gauge <a^+ a> for a given eigenvector (`eig_v`) in the RAD representation
/// calculated for the parameters in `prm`
fn ave_photon_pa(prm: &Parameters, eig_v:Array2<c64>) -> Array1<f64> {
    let i_m = iden(prm.n_kappa);

    let a = get_a(prm.nf, prm);

    let n_pa = kron(&i_m, &(&a.t()).dot(&a));
    let n_photons:Array2<c64> = (&eig_v).t().dot(&(&n_pa).dot(&eig_v));
    let n_ph_diag:Array1<f64> = n_photons.into_diag().iter().map(|&x| x.abs()).collect();

    n_ph_diag
}