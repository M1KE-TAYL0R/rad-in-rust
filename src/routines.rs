use csv::{ReaderBuilder,WriterBuilder};
use ndarray::*;
use std::{fs::File, result::Result, path::Path};
use ndarray_csv::{Array2Writer, Array2Reader, ReadError};

use crate::{parameters::*, graph_disp::*, solve_hamiltonian::*};


pub fn get_absorption(mut prm: Parameters, args: &Vec<String>) {

    for g_wc in prm.g_wc_grid.clone(){ 

        println!("Coupling = {}", g_wc);

        // Set coupling strength in prm
        prm.g_wc = g_wc;

        let mut data_total: Array3<f64> = Array3::zeros((prm.nk,prm.nk, prm.nf * prm.n_kappa));
        let mut data_total_c: Array3<f64> = Array3::zeros((prm.nk,prm.nk, prm.nf * prm.n_kappa));

        let k_ph_array = prm.k_points.clone();

        for k_ph in k_ph_array.iter().enumerate(){

            let wc_ph = (prm.wc_norm.powi(2) + (k_ph.1).powi(2)).sqrt();

            let mut data: Array2<f64> = Array2::zeros((prm.nk, prm.nf * prm.n_kappa + 1));
            let mut data_color: Array2<f64> = Array2::zeros((prm.nk, prm.nf * prm.n_kappa));

            (data, data_color) = rayon_dispatch(data, data_color, args, &prm.k_points, g_wc, Some(wc_ph));


            // let data_fname = filename(&prm, ".csv");
            // write_file(&data, &data_fname);

            // let color_fname = filename(&prm, "csv");
            // write_file(&data_color, &color_fname);
            let input = data.slice(s![..,1..]);
            data_total.slice_mut(s![k_ph.0,..,..]).assign(&input);
            data_total_c.slice_mut(s![k_ph.0,..,..]).assign(&data_color);

        }

        

    }
}




/// Generates a dispersion plot for each coupling strength in `prm.g_wc_grid`.
/// 
/// Iterates serially through each coupling strength, g_wc, and has the option to parallelize
/// the iteration over each k-point.
pub fn get_disps(mut prm: Parameters, args:&Vec<String> ) {
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
                (data, data_color) = rayon_dispatch(data, data_color, args, &prm.k_points, g_wc, None);

            }
            else {

                // Serial implementation of the iterator over k
                (data, data_color) = basic_dispatch(data, data_color, args, &prm.k_points, g_wc);

            }

            // Save the `data` and `data_color` arrays as .csv's
            write_file(&data, &data_fname);
            write_file(&data_color, &color_fname);
        }

        // Case where we are loading the data files
        else {
            // Load `read_data` from it's .csv
            let read_data: Result<Array2<f64>, ReadError> = 
                read_file(&data_fname,(prm.nk, prm.nf * prm.n_kappa + 1));

            // Save loaded data to `data` or print the reading error
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

        // Initialize image file name and plot data
        let image_fname = filename(&prm, "png");
        plot_data(&data, 30, &prm,&image_fname);
    }
}

/// Writes `data` to a .csv file of the name `fname`
pub fn write_file(data: &Array2<f64>, fname: &String) {
        
    let file = File::create(fname).unwrap();
    let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);
    writer.serialize_array2(&data).unwrap();

}

/// Reads a csv of the name `fname` and casts it to an `Array2` of shape `shape`
pub fn read_file(fname: &String, shape: (usize,usize)) -> Result<Array2<f64>, ReadError>{
    let file = File::open(fname).unwrap();
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
    let array_read: Result<Array2<f64>, ReadError> = reader.deserialize_array2(shape);

    array_read
}

/// Returns a unique file name based on the parameters in `prm` with file extension from `ext`
pub fn filename(prm: &Parameters, ext: &str) -> String{
    // let filename = format!("data/E_RAD_k{0:.3}_{1}_{2}_gwc{3:.7}_wc{4:.4}.dat",prm.k,prm.nf,prm.n_kappa,prm.g_wc,prm.wc_norm);

    let mut filename = format!("data/E_RAD_nk{0}_nf{1}_nkappa{2}_gwc{3:.7}_wc{4:.4}_kshift{5:.2}.",prm.nk,prm.nf,prm.n_kappa,prm.g_wc,prm.wc_norm,prm.k_shift);

    filename.push_str(ext);

    filename
}