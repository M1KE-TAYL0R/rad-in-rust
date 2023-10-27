#!/bin/bash 
#SBATCH -p action -A action 
#SBATCH -J 0_disp2d 
#SBATCH -o disp2d_output0.slurm 
#SBATCH -t 1-00:00:00 
#SBATCH --ntasks=24 
#SBATCH --mem-per-cpu=2GB 
 
/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust disp2d 1.0 -0.60205999132 -1.0 1 480 21 4