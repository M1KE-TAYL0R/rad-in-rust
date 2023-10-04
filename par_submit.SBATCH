#!/bin/bash 
#SBATCH -p action -A action 
#SBATCH -J 0_disp 
#SBATCH -o disp_output0.slurm 
#SBATCH -t 1-00:00:00 
#SBATCH -N 1 
#SBATCH --ntasks-per-node=24 
#SBATCH --mem 14GB 
 
/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust disp 1.0 -0.60205999132 -0.0 1 1200 21 7