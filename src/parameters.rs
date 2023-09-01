use std::{f64::consts::PI, usize};
use ndarray::*;

pub fn get_parameters() -> Parameters{
    let mut prm = Parameters{
        wc_norm: 1.0,
        wc: 1.0,
        ng: 128,
        g_wc: 0.0,
        g_wc_grid: array![0.0,0.1,0.2],
        omega: 1.0,
        xi_g: 0.0,
        nf: 5,
        n_cpus: 48,
        nk: 100,
        n_kappa: 51,
        n_kappa2: 11,
        k: 0.0,
        a_0: 4.0,
        z: 0.1278,
        r_0: 10.0,
        k_shift: 0.0,
        k_points: Array1::zeros(0),
        kappa_grid: Array1::zeros(0),
        kappa_grid2: Array1::zeros(0),
        m: 1.0,
        hbar: 1.0,
        load_existing: false
    };

    prm.k_points = Array1::linspace(-prm.a_0/PI + prm.k_shift, prm.a_0/PI + prm.k_shift, prm.nk);
    prm.kappa_grid  = PI / prm.a_0 *  Array1::linspace(prm.n_kappa  as f64 - 1.0, -(prm.n_kappa  as f64 - 1.0), prm.n_kappa );
    prm.kappa_grid2 = PI / prm.a_0 *  Array1::linspace(prm.n_kappa2 as f64 - 1.0, -(prm.n_kappa2 as f64 - 1.0), prm.n_kappa2);

    prm
}


pub struct Parameters {
    pub wc_norm: f64,
    pub wc: f64,
    pub g_wc: f64,
    pub ng: usize,
    pub g_wc_grid: Array1<f64>,
    pub nf: usize,
    pub n_cpus: usize,
    pub nk: usize,
    pub n_kappa: usize,
    pub n_kappa2: usize,
    pub k: f64,
    pub a_0: f64,
    pub z: f64,
    pub r_0: f64,
    pub k_shift: f64,
    pub k_points: Array1<f64>,
    pub kappa_grid: Array1<f64>,
    pub kappa_grid2: Array1<f64>,
    pub omega: f64,
    pub xi_g: f64,
    pub m: f64,
    pub hbar: f64,
    pub load_existing: bool
}

