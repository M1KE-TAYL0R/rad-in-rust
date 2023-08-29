use ndarray::*;
use ndarray_linalg::*;

fn main() {
    // let param = Parameters{omega:1.0, wc: 1.0};
    let b = get_b(3);
    // let a: = get_a(nf, param);
    println!("{}",b);
}

fn get_b(nf:usize) -> Array2<c64>{
    let mut b: Array2<c64> = Array2::zeros([nf,nf]);

    for m in 0..nf - 1 {
        let val = ((m + 1) as f64).sqrt();
        b[[m , m+1]] =c64::from(val);
    }

    b
}

fn get_a(nf:usize, param:&Parameters) -> Array2<c64>{
    let wc = param.wc;
    let omega = param.omega;

    let s_wc = (omega / wc).sqrt();

    let u = c64::from(0.5 * ( s_wc - 1.0 / s_wc ));
    let v = c64::from(0.5 * ( s_wc + 1.0 / s_wc ));

    let b = get_b(nf);

    u * &b - v * &b.t()
}

pub struct Parameters {
    wc_norm: f64,
    wc: f64,
    g_wc: f64,
    ng: usize,
    g_wc_grid: Array1<f64>,
    nf: usize,
    NCPUS: usize,
    nk: usize,
    n_kappa: usize,
    k: f64,
    a_0: f64,
    z: f64,
    r_0: f64,
    k_points: Array1<f64>,
    kappa_grid: Array1<f64>,
    kappa_grid2: Array1<f64>,
    omega: f64,
    xi_g: f64,
    m_0: f64,
    hbar: f64,
    load_existing: bool
}