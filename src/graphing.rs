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

pub fn plot_absorb(histogram:&Array2<f64>,prm: &Parameters, n_bins: usize, fname:&String) {
    let mut fig = Figure::new();
    fig.set_terminal("pngcairo size 1440,1080", fname)
    .set_pre_commands("set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#000000' behind");
    // fig.set_terminal("svg enhanced size 1440,1080", fname);

    let k_min = -PI / prm.a_0;
    let k_max = PI / prm.a_0;

    let mut x_min = k_min;
    if prm.near_edge{
        x_min = -0.4;
    }

    fig.axes2d()
    .set_x_range(AutoOption::Fix(x_min), AutoOption::Fix(-x_min))
    .set_y_range(AutoOption::Fix(0.0), AutoOption::Fix(prm.max_energy))
    .set_x_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    .set_y_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    .set_cb_ticks(Some((Auto,5)),&[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    .set_cb_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    .set_x_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    .set_y_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    .set_margins(&[MarginLeft(0.07),MarginBottom(0.07), MarginRight(0.87)])
    .set_border(true, &[Bottom,Right,Left,Top], &[Color("white")])
    .image(histogram, n_bins, n_bins, Some((k_min,0.0,k_max,prm.max_energy)), &[])
    ;

    let message = fig.save_to_png(fname, 1440, 1080);
    // let message_svg = fig.save_to_svg(fname, 1440, 1080);
    println!("{:?}", message);
}