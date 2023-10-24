use std::f64::consts::PI;
use ndarray::{*,linalg::kron};
use ndarray_linalg::{c64,*};
use statrs::function::gamma::gamma_ui;

use crate::parameters::*;

pub fn construct_pa_h(prm:&Parameters) -> Array2<c64> {
    // Define identities: 
    let i_m = iden(prm.n_kappa);
    let i_ph = iden(prm.nf);


    let k_e = get_t_pa(&prm);
    let v = get_v(&prm);
    let h_ph = get_h_ph(&prm);

    let mut h_total: Array2<c64> = Array2::zeros((prm.nf*prm.n_kappa,prm.nf*prm.n_kappa));

    h_total = h_total + kron(&i_m, &h_ph);
    h_total = h_total + k_e;
    h_total = h_total + kron(&v, &i_ph);

    h_total
}

fn get_t_pa(prm:&Parameters) -> Array2<c64>{
    
    let i_m = iden(prm.n_kappa);
    let i_ph = iden(prm.nf);

    let a = get_a_pa(prm.nf);
    let a_k = c64::from( prm.g_wc * (prm.hbar * prm.m * prm.wc).sqrt() );
    let vec_pot = a_k * (&a + &a.t());

    let p = to_c_2(Array2::from_diag(&(prm.kappa_grid.clone() + prm.k)));

    let mut p_new = c64::from(prm.hbar) * kron(&p, &i_ph);
    p_new = p_new - kron(&i_m,  &(a.t().dot(&a))) * c64::from(prm.hbar*(prm.k - prm.k_shift));
    p_new = p_new - kron(&i_m, &vec_pot);

    let k_e = &p_new.dot(&p_new) / (2.0 * prm.m);

    k_e
}

/// Calculates V_RAD for a given set of parameters `prm`
fn get_v(prm:&Parameters) -> Array2<c64>{
    // Generate V
    let mut v: Array2<c64> = Array2::zeros((prm.n_kappa,prm.n_kappa));

    let kappa_max = prm.kappa_grid.norm_max();
    
    for (l, k2) in prm.kappa_grid2.iter().enumerate() {
        for (m, k1) in prm.kappa_grid.iter().enumerate(){
            if ((k1 + k2).abs() <= kappa_max) && !(k2.abs() <= f64::EPSILON) {
                let shift_ind = (prm.n_kappa2 - 1) / 2;
                v[[m,m + l - shift_ind]] =  c64::from(- prm.z / 2.0 / PI * gamma_ui(1e-10, (k2.re() / 2.0 / prm.r_0).powi(2)));
            }
        }
    }
    v
}

/// Returns H_ph for a given set of parameters `prm`
fn get_h_ph(prm:&Parameters) -> Array2<c64>{

    let range = to_c_1(Array1::linspace(0.0, (prm.nf - 1) as f64 , prm.nf));

    let h_ph: Array2<c64> = Array2::from_diag(&range) * prm.wc;

    h_ph
}

/// Calculates the Bogoliubov transformed b operator matrix
pub fn get_a_pa(nf:usize) -> Array2<c64>{
    let mut a: Array2<c64> = Array2::zeros([nf,nf]);

    for m in 0..nf - 1 {
        let val = ((m + 1) as f64).sqrt();
        a[[m , m+1]] =c64::from(val);
    }

    a
}

