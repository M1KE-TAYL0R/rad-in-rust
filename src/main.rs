use csv::*;
use ndarray::*;
use std::fs::File;
use ndarray_csv::Array2Writer;
use rayon::{iter::IntoParallelRefMutIterator, prelude::IndexedParallelIterator};
use rayon::iter::ParallelIterator;
use std::time::Instant;
use std::env;

mod build_hamiltonian;
use build_hamiltonian::*;

pub mod parameters;
use parameters::*;

mod solve_hamiltonian;
use solve_hamiltonian::*;

mod graph_disp;
use graph_disp::*;

fn main() {
    let now = Instant::now();
    let mut prm = get_parameters();

    // Read command line arguments: log_g_min, log_g_max, ng, nk, n_kappa, nf
    let args: Vec<String> = env::args().collect();

    if args.len() == 6 {
        let log_g_min = args[0].parse::<f64>().unwrap();
        let log_g_max = args[1].parse::<f64>().unwrap();
        prm.ng = args[2].parse::<usize>().unwrap();
        prm.g_wc_grid = Array1::linspace(log_g_min, log_g_max, prm.ng).map(|x| libm::exp10(*x as f64));

        prm.nk = args[3].parse::<usize>().unwrap();
        prm.n_kappa = args[4].parse::<usize>().unwrap();
        prm.nf = args[5].parse::<usize>().unwrap();
    }
    else{
        panic!();
    }

    let k_points = prm.k_points.clone();
    let g_wc_grid = prm.g_wc_grid.clone();
 
    // println!("{}",prm.kappa_grid2);

    for g_wc in &g_wc_grid{
        println!("Coupling = {}", g_wc);

        let rayon = true;

        let mut data: Array2<f64> = Array2::zeros((prm.nk, prm.nf * prm.n_kappa + 1));

        if rayon {
            ////////////////// Rayon Version //////////////////////////

            data.slice_mut(s![..,0]).assign(&k_points);

            prm.g_wc = *g_wc;

            let mut data_vec: Vec<Array1<f64>> = Vec::new();

            for ijk in 0 .. prm.nk {
                data_vec.push(data.slice(s![ijk,1..]).to_owned());
            }

            data_vec.par_iter_mut().enumerate().for_each( |(k, col)| {
                
                let mut prm_k = get_parameters();

                prm_k.k = k_points[k];

                prm_k.wc = (prm_k.wc_norm.powi(2) + (k_points[k] - prm_k.k_shift).powi(2)).sqrt();

                (prm_k.omega, prm_k.xi_g) = get_couplings(&prm_k);

                (*col,_) = solve_h(&prm_k);
            
            });

            for ijk in 0 .. prm.nk {
                data.slice_mut(s![ijk,1..]).assign( &data_vec[ijk] );
            }
        }
        else {
            /////////////// Basic Version /////////////////////////////
            data.slice_mut(s![..,0]).assign(&k_points);
            
            for k in k_points.iter().enumerate(){
                prm.k = *k.1;
                prm.g_wc = *g_wc;
                prm.wc = (prm.wc_norm.powi(2) + (k.1 - prm.k_shift).powi(2)).sqrt();

                (prm.omega, prm.xi_g) = get_couplings(&prm);

                let (eig_e,_) = solve_h(&prm);
                data.slice_mut(s![k.0,1..]).assign(&eig_e);
            }
        }


        let data_fname = filename(&prm, "csv");
        let message = write_file(&data, &data_fname);
        println!("{:?}", message);

        let image_fname = filename(&prm, "png");
        plot_data(&data, 30, &prm,&image_fname);
    }

    println!("Time elapsed = {} sec", now.elapsed().as_secs());
}


fn write_file(data: &Array2<f64>, fname: &String) -> Result<()> {
    
    {
        let file = File::create(fname)?;
        let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);
        writer.serialize_array2(&data)?;
    }

    Ok(())
}
