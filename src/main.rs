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

/*

Rust implementation of the RAD Hamiltonian beyong the long-wavelength approximation
for k modes coupled to a single particle. This program solves the dispersion 
relations of a single electron in an infinite periodic erf() potential for multiple
coupling strengths and k-points with the option of displacing the photonic k from 
the matter k using the k-shift parameter.

To run perform the following commands:

cargo build --release
./target/release/rad-in-rust {log_g_min} {log_g_max} {ng} {nk} {n_kappa} {nf}

For use on bluehive, the parallel.py script dispatches this routine accross multiple
nodes.

*/


/// Main function which iterates serially through each coupling strength, g_wc, and
/// has the option to parallelize the iteration over each k-point.
fn main() {
    
    // Start the timer to determine
    let now = Instant::now();

    // Read command line arguments: log_g_min, log_g_max, ng, nk, n_kappa, nf
    let args: Vec<String> = env::args().collect();
    
    // Initialize the prm struct based off the command-line arguments
    let mut prm = get_parameters(&args);

    // Loop over all coupling strengths in the array prm.g_wc_grid
    for g_wc in prm.g_wc_grid.clone(){

        println!("Coupling = {}", g_wc);

        // Set coupling strength in prm
        prm.g_wc = g_wc;

        // Flag to turn on/off the rayon parallelization implementation
        let rayon: bool = true;

        // Initialize the arrays to store the energies (data) and photon numbers (data_color)
        // in the same scope as the saving function write_file()
        let mut data: Array2<f64> = Array2::zeros((prm.nk, prm.nf * prm.n_kappa + 1));
        let mut data_color: Array2<f64> = Array2::zeros((prm.nk, prm.nf * prm.n_kappa));

        // Initialize the data filenames in the same scope as the saving function write_file()
        let data_fname = filename(&prm, "csv");
        let color_fname: String = filename(&prm, "_color.csv");

        // Determines if we need to calculate the eigenenergies or if they were previously 
        // calculated and can/should be loaded
        if !(prm.load_existing && Path::new(&data_fname).exists() && Path::new(&color_fname).exists()){

            // Decides to either do the parallel or serial implementation based on the bool, rayon
            if rayon {

                // Rayon (parallel) implementation of the iterator over k
                (data, data_color) = rayon_dispatch(data, data_color, &args, &prm.k_points, g_wc);

            }
            else {

                // Serial implementation of the iterator over k
                (data, data_color) = basic_dispatch(data, data_color, &args, &prm.k_points, g_wc);

            }


            write_file(&data, &data_fname);
            write_file(&data_color, &color_fname);
        }

        else {
            let read_data: Result<Array2<f64>, ReadError> = 
                read_file(&data_fname,(prm.nk, prm.nf * prm.n_kappa + 1));

            match read_data {
                Ok(v) => data = v,
                Err(e) => println!("Reading error: {e:?}"),
            }

            // let read_color: Result<Array2<f64>, ReadError> = 
            //     read_file(&color_fname,(prm.nk, prm.nf * prm.n_kappa ));

            // match read_color {
            //     Ok(v) => data_color = v,
            //     Err(e) => println!("Reading error: {e:?}"),
            // }
        }

        let image_fname = filename(&prm, "png");
        plot_data(&data, 30, &prm,&image_fname);
    }

    println!("Time elapsed = {} sec", now.elapsed().as_secs());
}

fn rayon_dispatch(mut data: Array2<f64>, mut data_color: Array2<f64>, args:&Vec<String>, k_points: &Array1<f64>, g_wc: f64) -> (Array2<f64>, Array2<f64>) {
    let prm = get_parameters(args);

    data.slice_mut(s![..,0]).assign(&k_points);

    let mut data_vec: Vec<(Array1<f64>,Array1<f64>)> = vec![(Array1::zeros(prm.nf * prm.n_kappa + 1),Array1::zeros(prm.nf * prm.n_kappa + 1)); prm.nk];

    data_vec.par_iter_mut().progress().enumerate().for_each( |(k, col)| {
        
        let mut prm_k = get_parameters(&args);

        prm_k.k = k_points[k];
        prm_k.g_wc = g_wc;

        prm_k.wc = (prm_k.wc_norm.powi(2) + (k_points[k] - prm_k.k_shift).powi(2)).sqrt();

        (prm_k.omega, prm_k.xi_g) = get_couplings(&prm_k);

        *col = solve_h(&prm_k);
    
    });

    for ijk in 0 .. prm.nk {
        data.slice_mut(s![ijk,1..]).assign( &(data_vec[ijk].0) );
        data_color.slice_mut(s![ijk,..]).assign( &(data_vec[ijk].1) );
    }

    (data, data_color)
}

fn basic_dispatch(mut data: Array2<f64>, mut data_color: Array2<f64>, args:&Vec<String>, k_points: &Array1<f64>, g_wc: f64) -> (Array2<f64>, Array2<f64>){

    let mut prm = get_parameters(args);

    prm.g_wc = g_wc;

    data.slice_mut(s![..,0]).assign(k_points);
                
    for k in k_points.iter().progress().enumerate(){
        prm.k = *k.1;
        prm.wc = (prm.wc_norm.powi(2) + (k.1 - prm.k_shift).powi(2)).sqrt();

        (prm.omega, prm.xi_g) = get_couplings(&prm);

        let (eig_e,n_pa) = solve_h(&prm);
        data.slice_mut(s![k.0,1..]).assign(&eig_e);
        data_color.slice_mut(s![k.0,1..]).assign(&n_pa);
    }

    (data, data_color)
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

/// Generates a unique file name based on the parameters in `prm` with file extension from `ext`
fn filename(prm: &Parameters, ext: &str) -> String{
    // let filename = format!("data/E_RAD_k{0:.3}_{1}_{2}_gwc{3:.7}_wc{4:.4}.dat",prm.k,prm.nf,prm.n_kappa,prm.g_wc,prm.wc_norm);

    let mut filename = format!("data/E_RAD_nk{0}_nf{1}_nkappa{2}_gwc{3:.7}_wc{4:.4}_kshift{5:.2}.",prm.nk,prm.nf,prm.n_kappa,prm.g_wc,prm.wc_norm,prm.k_shift);

    filename.push_str(ext);

    filename
}