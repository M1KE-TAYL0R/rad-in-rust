#!/bin/bash 
#SBATCH -p standard 
#SBATCH -J 7_absorb 
#SBATCH -o absorb_output7.slurm 
#SBATCH -t 1-00:00:00 
#SBATCH -N 1 
#SBATCH --ntasks-per-node=24 
#SBATCH --mem 60GB 
 
cargo build --release 
 
/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust absorb 0.1 -0.25 -0.0 1 360 31 10