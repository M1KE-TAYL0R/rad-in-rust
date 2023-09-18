#!/bin/bash 
#SBATCH -p standard 
#SBATCH -J 9_absorb 
#SBATCH -o absorb_output9.slurm 
#SBATCH -t 1-00:00:00 
#SBATCH -N 1 
#SBATCH --ntasks-per-node=36 
#SBATCH --mem 60GB 
 
cargo build --release 
 
/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust absorb 0.08 -0.6500000000000001 -0.5 1 360 101 7