#!/bin/bash 
#SBATCH -p preempt 
#SBATCH -J 0_rad_rust 
#SBATCH -o output0.slurm 
#SBATCH -t 2-00:00:00 
#SBATCH -N 1 
#SBATCH --ntasks-per-node=24 
#SBATCH --mem 60GB 
 
cargo build --release 
 
/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust -1.0 1.0 4 120 31 2