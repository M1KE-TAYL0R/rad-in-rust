#!/bin/bash
#SBATCH -p action -A action
#SBATCH -J rad_rust
#SBATCH -o output.slurm
#SBATCH -t 2-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mem 60GB

cargo build --release
target/release/rad-in-rust -1.0 0.0 3 5 31 5