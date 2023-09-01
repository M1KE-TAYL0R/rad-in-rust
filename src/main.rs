use csv::*;
use ndarray::*;
use std::fs::File;
use ndarray_csv::Array2Writer;

mod build_hamiltonian;
use build_hamiltonian::*;

pub mod parameters;
use parameters::*;

mod solve_hamiltonian;
use solve_hamiltonian::*;

mod graph_disp;
use graph_disp::*;

fn main() {
    let mut prm = get_parameters();

    let k_points = prm.k_points.clone();
    let g_wc_grid = prm.g_wc_grid.clone();
 
    // println!("{}",prm.kappa_grid2);

    for g_wc in &g_wc_grid{
        println!("Coupling = {}", g_wc);
        let mut data: Array2<f64> = Array2::zeros((prm.nk, prm.nf * prm.n_kappa + 1));
        data.slice_mut(s![..,0]).assign(&k_points);
        
        for k in k_points.iter().enumerate(){
            prm.k = *k.1;
            prm.g_wc = *g_wc;
            prm.wc = (prm.wc_norm.powi(2) + (k.1 - prm.k_shift).powi(2)).sqrt();

            (prm.omega, prm.xi_g) = get_couplings(&prm);

            let (eig_e,_) = solve_h(&prm);
            data.slice_mut(s![k.0,1..]).assign(&eig_e);
        }

        let fname = filename(&prm);
        _ = write_file(&data, &fname);
        plot_data(&data, 5, &prm,&fname);
    }


}


fn write_file(data: &Array2<f64>, fname: &String) -> Result<()> {
    
    {
        let file = File::create(fname)?;
        let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);
        writer.serialize_array2(&data)?;
    }

    Ok(())
}
