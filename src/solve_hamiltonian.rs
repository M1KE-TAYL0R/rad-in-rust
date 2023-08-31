use ndarray::{*,linalg::kron};
use ndarray_linalg::{c64,*};

use crate::parameters::*;

pub fn filename(prm: &Parameters) -> String{
    let filename = format!("data/E_RAD_k{0:.3}_{1}_{2}_gwc{3:.7}_wc{4:.4}.dat",prm.k,prm.nf,prm.n_kappa,prm.g_wc,prm.wc_norm);
    filename
}

pub fn solve_h (mut prm: Parameters, k: f64, g_wc: f64) {
    prm.k = k;
    prm.g_wc = g_wc;
    prm.wc = (prm.wc_norm.powi(2) + (k - prm.k_shift).powi(2)).sqrt();

    let filename = filename(&prm);

    

}