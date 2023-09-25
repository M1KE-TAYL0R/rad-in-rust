use std::{usize,f64::consts::PI};
use ndarray::*;

/// From the command line arguments `args` this populates a `Parameters` struct and returns it
pub fn get_parameters(args: &Vec<String>) -> Parameters{
    
    let mut prm = Parameters{
        wc_norm: 0.1,
        wc: 0.0,
        ng: 128,
        g_wc: 0.0,
        g_wc_grid: Array1::zeros(0),
        omega: 1.0,
        xi_g: 0.0,
        nf: 5,
        n_cpus: 48,
        nk: 240,
        n_kappa: 64,
        n_kappa2: 11,
        k: 0.0,
        k_ph: 0.0,
        a_0: 4.0,
        z: 0.1278,
        r_0: 10.0,
        k_shift: 0.0,
        k_points: Array1::zeros(0),
        kappa_grid: Array1::zeros(0),
        kappa_grid2: Array1::zeros(0),
        m: 1.0,
        hbar: 1.0,
        load_existing: false,
        max_energy: 1.5,
        k_ph_factor: 15,
        near_edge: true
    };

    // prm.k_points = Array1::linspace(-prm.a_0/PI + prm.k_shift, prm.a_0/PI + prm.k_shift, prm.nk);
    // prm.kappa_grid  = PI / prm.a_0 *  Array1::linspace(prm.n_kappa  as f64 - 1.0, -(prm.n_kappa  as f64 - 1.0), prm.n_kappa );
    // prm.kappa_grid2 = PI / prm.a_0 *  Array1::linspace(prm.n_kappa2 as f64 - 1.0, -(prm.n_kappa2 as f64 - 1.0), prm.n_kappa2);

    // prm.g_wc_grid = Array1::linspace(-1.0, 2.0, prm.ng).map(|x| libm::exp10(*x as f64));  

    // println!("{}", args[1]);

    prm.wc_norm = args[1].parse::<f64>().unwrap();

    let log_g_min = args[2].parse::<f64>().unwrap();
    let log_g_max = args[3].parse::<f64>().unwrap();
    prm.ng = args[4].parse::<usize>().unwrap();
    let temp = Array1::linspace(log_g_min, log_g_max, prm.ng + 1).map(|x| libm::exp10(*x as f64));
    prm.g_wc_grid = temp.slice(s![..-1]).to_owned();

    prm.nk = args[5].parse::<usize>().unwrap();
    prm.n_kappa = args[6].parse::<usize>().unwrap();
    prm.nf = args[7].parse::<usize>().unwrap();

    // prm.k_shift = -PI/prm.a_0; // Added for debugging!

    if prm.near_edge{
        prm.k_points = Array1::linspace(PI/prm.a_0 - 0.03, PI/prm.a_0 + 0.03, prm.nk); // k-points near the K-point
    }
    else{
        prm.k_points = Array1::linspace(-PI/prm.a_0 + prm.k_shift, PI/prm.a_0 + prm.k_shift, prm.nk);
    }

    prm.kappa_grid  = PI / prm.a_0 *  Array1::linspace(prm.n_kappa  as f64 - 1.0, -(prm.n_kappa  as f64 - 1.0), prm.n_kappa );
    prm.kappa_grid2 = PI / prm.a_0 *  Array1::linspace(prm.n_kappa2 as f64 - 1.0, -(prm.n_kappa2 as f64 - 1.0), prm.n_kappa2);

    prm
}

/// Struct containing all of the necessary infomation to solve a dispersion plot
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
    pub k_ph: f64,
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
    pub load_existing: bool,
    pub max_energy: f64,
    pub k_ph_factor: usize,
    pub near_edge: bool
}

