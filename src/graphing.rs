use std::f64::consts::PI;
use gnuplot::*;
use ndarray::{Array2,Array1,s,Array3};
use plotters::{prelude::*, style::Color, coord::{types::RangedCoordf64, cartesian::Cartesian3d}};
use statrs::statistics::Statistics;

use crate::parameters::*;

/// Uses the `gnuplot` crate to plot the dispersion plots.
/// Currently cannot graph the photon number yet.
pub fn _plot_disp(data:&Array2<f64>, n_states:usize, prm: &Parameters,fname:&String) {
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

pub fn plotters_disp(data:&Array2<f64>,data_c:&Array2<f64>, n_states:usize, prm: &Parameters,fname:&str) -> Result<(), Box<dyn std::error::Error>>{
    let x_max = (PI / prm.a_0) as f32;
    let y_max = 5.0 as f32;
    let c_max = data_c.max() as f32;

    let scale_factor = 1;

    let zpe = data.column(1).to_vec().into_iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let zpe_vec: Array1<f64> = Array1::ones(prm.nk)*zpe;

    let x = data.slice(s![..,0]);

    let root = SVGBackend::new(fname, (1440*scale_factor,1080*scale_factor)).into_drawing_area();

    let original_style = ShapeStyle {
        color: BLACK.into(),
        filled: false,
        stroke_width: 2*scale_factor,
    };

    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
        .x_label_area_size(70*scale_factor)
        .y_label_area_size(100*scale_factor)
        .build_cartesian_2d(-x_max..x_max, 0.0 .. y_max)?;

    chart.configure_mesh()
    .x_label_style(("helvetica", 50*scale_factor))
    .y_label_style(("helvetica", 50*scale_factor))
    .disable_mesh()
    .set_all_tick_mark_size(20*scale_factor)
    .bold_line_style(original_style)
    .x_label_formatter(&|v| format!("{:.1}", v))
    .y_label_formatter(&|v| format!("{:.1}", v))
    .draw()?;


    for m in 0 .. n_states {
        let col = data.column(m+1);
        let color = data_c.column(m);
        let y = col.to_owned() - &zpe_vec;

        // chart.draw_series(x.iter().enumerate().map(|(ind, el)| {
        //     if ind != 0 && y[ind] < (y_max as f64){ 
        //         let temp = [(*el as f32, y[ind] as f32), (x[ind -1] as f32, y[ind-1] as f32)];
        //         let temp_c = color[ind];
        //         let c = VulcanoHSL::get_color_normalized(temp_c as f32, 0.0, c_max);
        //         PathElement::new(temp, c.stroke_width(4))
        //     }
        //     else {
        //         PathElement::new([(0.0 as f32,0.0 as f32), (0.0 as f32,0.0 as f32)], WHITE)
        //     }
        // }))?;

        chart.draw_series(x.iter().enumerate().map(|(ind, el)| {
            if y[ind] < (y_max as f64){
            let temp = (*el as f32, y[ind] as f32);
            EmptyElement::at(temp)
                + Circle::new((0,0), 3*scale_factor, VulcanoHSL::get_color_normalized(color[ind] as f32, 0.0, c_max).filled())
            }
            else {
                EmptyElement::at((0.,0.)) + Circle::new((0,0), 1, VulcanoHSL::get_color(0.0))
            }
        }))?;

    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", fname);

    Ok(())
}

pub fn plotters_disp_2d(data:&Array3<f64>,data_c:&Array3<f64>,n_states:usize, prm: &Parameters,fname:&str) -> Result<(), Box<dyn std::error::Error>>{
    let x_max = PI / prm.a_0;
    let z_max = PI / prm.a_0;
    let y_max = 2.0;
    let _c_max = data_c.max();

    // let x_array = Array1::linspace(-x_max, x_max, prm.nk * 2);
    // let z_array = Array1::linspace(-x_max, x_max, prm.nk * 2);

    let scale_factor = 2;

    // let zpe = data.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    let gif = false;

    if gif {
        let fps = 30;
        let length = 4;
        let area = BitMapBackend::gif(
            "disp/animated.gif", 
            (1080*scale_factor, 1080*scale_factor), 
            1000/fps
        )?.into_drawing_area();
            
        for i in 0..=fps*length {
            area.fill(&WHITE)?;
            
            let mut chart = ChartBuilder::on(&area)
                .margin(20)
                // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
                .x_label_area_size(70*scale_factor)
                .y_label_area_size(100*scale_factor)
                .build_cartesian_3d(-x_max..x_max, 0.0 .. y_max, - z_max .. z_max)?;

            chart.with_projection(|mut pb| {
                pb.pitch = 0.2;
                pb.yaw = (i as f64) / ((fps * length) as f64) * PI ;
                pb.scale = 0.8;
                pb.into_matrix()
            });

            chart
            // .set_3d_pixel_range((1440,1440,1440))
            .configure_axes()
            .label_style(("helvetica", 25*scale_factor))
            .draw()?;

            draw_s_series(chart, n_states, data, data_c, prm, (x_max,y_max,z_max));

            area.present().unwrap();
        }

        Ok(())
    }

    else {
        // let root = SVGBackend::new(fname, (1440*scale_factor,1080*scale_factor)).into_drawing_area();
        let root = BitMapBackend::new(fname, (1440*scale_factor,1080*scale_factor)).into_drawing_area();

        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(&root)
            .margin(20)
            // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
            .x_label_area_size(70*scale_factor)
            .y_label_area_size(100*scale_factor)
            .build_cartesian_3d(-x_max..x_max, 0.0 .. y_max, - z_max .. z_max)?;

        chart.with_projection(|mut pb| {
            pb.pitch = 0.2;
            pb.yaw = 2.2;
            pb.scale = 0.95;
            pb.into_matrix()
        });
        // chart.with_projection(|mut pb| {
        //     pb.pitch = 0.0;
        //     pb.yaw = PI / 2.0;
        //     pb.scale = 0.95;
        //     pb.into_matrix()
        // });

        chart
        // .set_3d_pixel_range((1440,1440,1440))
        .configure_axes()
        .label_style(("helvetica", 25*scale_factor))
        .draw()?;

        draw_s_series(chart, n_states, data, data_c, prm, (x_max,y_max,z_max));

        root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
        println!("Result has been saved to {}", fname);

        Ok(())
    }
}

fn draw_s_series(mut chart:ChartContext<'_, BitMapBackend<'_>, Cartesian3d<RangedCoordf64, RangedCoordf64, RangedCoordf64>>, n_states:usize, data:&Array3<f64>, data_c:&Array3<f64>, prm: &Parameters, maxes: (f64,f64,f64)) {
    let x_max = maxes.0;
    let z_max = maxes.2;

    let x_array = Array1::linspace(-x_max, x_max, prm.nk + 1);
    let z_array = Array1::linspace(-x_max, x_max, prm.nk + 1);

    let zpe = data.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    for m in 0 .. n_states {
        chart.draw_series(SurfaceSeries::xoz(
            x_array.iter().map(|el| *el),
            z_array.iter().map(|el| *el),
            |x:f64,z:f64| {
                let i = ((x + x_max) / x_max / 2.0 * (prm.nk - 1) as f64) as usize;
                let j = ((z + z_max) / z_max / 2.0 * (prm.nk - 1) as f64) as usize;
                // println!("i {} j {} m{} data{}", i, j, m, data_c[[i,j,m]]);
                data[[i,j,m]] - *zpe
            } ).style_func(
                &|y: &f64| 
                {
                    let _val = *y;
                    if *y > maxes.1{
                        HSLColor(0.0, 0.0, 0.0).mix(0.0).filled()
                    }
                    else {
                        let col = m as f64 / n_states as f64;
                        // println!("{}", col);
                        // HSLColor(0.333, col as f64, 0.7).mix(0.5).filled()
                        HSLColor(0.6666, col as f64, 0.666).mix(0.5).filled()
                    }
                }
            )
        ).unwrap();

        // chart.draw_series(LineSeries::new(
        //     (0..prm.nk).map(|ind| (x_array[0],data[[0,ind,m]],x_array[ind])),
        //     (&BLACK).stroke_width(3)
        // )).unwrap();

        chart.draw_series( (0..prm.nk).map({ |ind|
            if (data[[ind,ind,m]] - zpe) < maxes.1 {
                EmptyElement::at((x_array[ind],data[[ind,ind,m]] - zpe,x_array[ind])) + Circle::new((0,0), 6, VulcanoHSL::get_color_normalized(data_c[[ind,ind,m]] as f32, 0.0, 1.0).filled())
            }
            else {
                EmptyElement::at((0.,0.,0.)) + Circle::new((0,0), 1, &WHITE.mix(0.0))
            }
        })).unwrap();
    }

    for m in 0 .. n_states {
        chart.draw_series( (0..prm.nk).map({ |ind|
            if (data[[0,ind,m]] - zpe) < maxes.1 {
                EmptyElement::at((-x_array[0],data[[0,ind,m]] - zpe,x_array[ind])) + Circle::new((0,0), 3, (&BLACK).filled())
            }
            else {
                EmptyElement::at((0.,0.,0.)) + Circle::new((0,0), 1, &WHITE.mix(0.0))
            }
        })).unwrap();

        chart.draw_series( (0..prm.nk).map({ |ind|
            if (data[[ind,0,m]] - zpe) < maxes.1 {
                EmptyElement::at((x_array[ind],data[[ind,0,m]] - zpe,x_array[0])) + Circle::new((0,0), 3, (&BLACK).filled())
            }
            else {
                EmptyElement::at((0.,0.,0.)) + Circle::new((0,0), 1, &WHITE.mix(0.0))
            }
        })).unwrap();
    }

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