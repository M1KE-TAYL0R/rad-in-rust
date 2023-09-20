#!/bin/bash 
#SBATCH -p action -A action 
#SBATCH -J 0_absorb 
#SBATCH -o absorb_outputK0.slurm 
#SBATCH -t 1-00:00:00 
#SBATCH -N 1 
#SBATCH --ntasks-per-node=24 
#SBATCH --mem 24GB 
 
cargo build --release 
 
/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust absorb 0.1 -1.0 -0.0 1 360 21 5