#!/bin/bash 
#SBATCH -p action -A action 
#SBATCH -J 7_absorb 
#SBATCH -o absorb_output07.slurm 
#SBATCH -t 1-00:00:00 
#SBATCH -N 1 
#SBATCH --ntasks-per-node=24 
#SBATCH --mem 14GB 
 
cargo build --release 
 
/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust absorb 0.15 -0.25 -0.0 2 2880 21 10