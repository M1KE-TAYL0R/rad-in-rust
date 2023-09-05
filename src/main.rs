use csv::*;
use ndarray::*;
use std::fs::File;
use ndarray_csv::{Array2Writer, Array2Reader, ReadError};
use rayon::{iter::IntoParallelRefMutIterator, prelude::IndexedParallelIterator};
use rayon::iter::ParallelIterator;
use std::time::Instant;
use std::env;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use std::path::Path;
use std::result::Result;

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

    // Read command line arguments: log_g_min, log_g_max, ng, nk, n_kappa, nf
    let args: Vec<String> = env::args().collect();
    
    let mut prm = get_parameters(&args);
  
    let k_points = prm.k_points.clone();
    let g_wc_grid = prm.g_wc_grid.clone();
 
    // println!("{}",prm.kappa_grid2);

    for g_wc in &g_wc_grid{
        println!("Coupling = {}", g_wc);

        let rayon = false;

        let mut data: Array2<f64> = Array2::zeros((prm.nk, prm.nf * prm.n_kappa + 1));

        let data_fname = filename(&prm, "csv");

        prm.g_wc = *g_wc;

        if !(prm.load_existing && Path::new(&data_fname).exists()){

            if rayon {
                ////////////////// Rayon Version //////////////////////////
            
                data.slice_mut(s![..,0]).assign(&k_points);

                let mut data_vec: Vec<Array1<f64>> = Vec::new();

                for ijk in 0 .. prm.nk {
                    data_vec.push(data.slice(s![ijk,1..]).to_owned());
                }

                data_vec.par_iter_mut().progress().enumerate().for_each( |(k, col)| {
                    
                    let mut prm_k = get_parameters(&args);

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
                
                for k in k_points.iter().progress().enumerate(){
                    prm.k = *k.1;
                    prm.wc = (prm.wc_norm.powi(2) + (k.1 - prm.k_shift).powi(2)).sqrt();

                    (prm.omega, prm.xi_g) = get_couplings(&prm);

                    let (eig_e,_) = solve_h(&prm);
                    data.slice_mut(s![k.0,1..]).assign(&eig_e);
                }
            }


            write_file(&data, &data_fname);
        }

        else {
            let read_data = read_file(&data_fname,(prm.nk, prm.nf * prm.n_kappa + 1));

            match read_data {
                Ok(v) => data = v,
                Err(e) => println!("Reading error: {e:?}"),
            }

            // data = read_file(&data_fname).unwrap();
        }

        let image_fname = filename(&prm, "png");
        plot_data(&data, 30, &prm,&image_fname);
    }

    println!("Time elapsed = {} sec", now.elapsed().as_secs());
}



fn write_file(data: &Array2<f64>, fname: &String) {
        
    let file = File::create(fname).unwrap();
    let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);
    writer.serialize_array2(&data).unwrap();

}

fn read_file(fname: &String, shape: (usize,usize)) -> Result<Array2<f64>, ReadError>{
    let file = File::open(fname).unwrap();
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
    let array_read: Result<Array2<f64>, ReadError> = reader.deserialize_array2(shape);

    array_read
}