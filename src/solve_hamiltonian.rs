use ndarray::{*,linalg::kron};
use ndarray_linalg::{c64,*};

use crate::{parameters::*, build_hamiltonian::*};

pub fn filename(prm: &Parameters, ext: &str) -> String{
    // let filename = format!("data/E_RAD_k{0:.3}_{1}_{2}_gwc{3:.7}_wc{4:.4}.dat",prm.k,prm.nf,prm.n_kappa,prm.g_wc,prm.wc_norm);

    let mut filename = format!("data/E_RAD_nk{0}_nf{1}_nkappa{2}_gwc{3:.7}_wc{4:.4}_kshift{5:.2}.",prm.nk,prm.nf,prm.n_kappa,prm.g_wc,prm.wc_norm,prm.k_shift);

    filename.push_str(ext);

    filename
}

pub fn solve_h (prm: &Parameters) -> (Array1<f64>, Array1<f64>){

    let h = construct_h_total(prm);

    let (eig_e, eig_v) = h.eigh(UPLO::Upper).unwrap();

    let n_pa_diag = ave_photon_pa(&prm, eig_v);

    (eig_e,n_pa_diag)
}

fn ave_photon_pa(prm: &Parameters, eig_v:Array2<c64>) -> Array1<f64> {
    let i_m = iden(prm.n_kappa);

    let a = get_a(prm.nf, prm);

    let n_pa = kron(&i_m, &(&a.t()).dot(&a));
    let n_photons:Array2<c64> = (&eig_v).t().dot(&(&n_pa).dot(&eig_v));
    let n_ph_diag:Array1<f64> = n_photons.into_diag().iter().map(|&x| x.abs()).collect();

    n_ph_diag
}