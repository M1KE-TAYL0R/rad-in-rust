use std::f64::consts::PI;
use gnuplot::*;
use ndarray::{Array2,Array1,s};

use crate::parameters::*;

/// Uses the `gnuplot` crate to plot the dispersion plots.
/// Currently cannot graph the photon number yet.
pub fn plot_disp(data:&Array2<f64>, n_states:usize, prm: &Parameters,fname:&String) {
    let mut fig = Figure::new();
    fig.set_terminal("pngcairo size 1440,1080", fname);

    let x_max = PI / prm.a_0;

    let zpe = data.column(1).to_vec().into_iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let zpe_vec: Array1<f64> = Array1::ones(prm.nk)*zpe;

    let x = data.slice(s![..,0]);

    for m in 0..n_states {
        let col = data.column(m+1);
        let y = col.to_owned() - &zpe_vec;

        fig.axes2d().lines(x,y, &[Color("black"),LineWidth(1.5)])
        .set_x_range(AutoOption::Fix(-x_max + prm.k_shift), AutoOption::Fix(x_max + prm.k_shift))
        .set_y_range(AutoOption::Fix(0.0), AutoOption::Fix(5.0));
    }

    let message = fig.save_to_png(fname, 1440, 1080);
    println!("{:?}", message);
}

pub fn plot_absorb(histogram:&Array2<f64>,nk: usize, n_e_bins: usize, fname:&String) {
    let mut fig = Figure::new();
    fig.set_terminal("pngcairo size 1440,1080", fname);

    fig.axes2d().image(histogram, n_e_bins, nk, None, &[]);

    let message = fig.save_to_png(fname, 1440, 1080);
    println!("{:?}", message);
}