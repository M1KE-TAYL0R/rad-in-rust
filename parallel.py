import os
import numpy as np

# Compile on head node
os.system("cargo build --release")

# Simulation parameters
n_sections = 1
g_min_log =  -0.60205999132
# g_min_log =  -15.0
g_max_log =  -1.0
log_g_bounds = np.linspace(g_min_log,g_max_log,n_sections+1)
ng = 1
nf = 4
nk = 480
n_kappa = 21
wc_norm = 1.0
routine = "disp2d"

print(f"Running rad-in-rust with {ng*n_sections} couplings and {nk} k-points")

# sbatch parameters
# node = "-p preempt"
node = "-p action -A action"
time = "1-00:00:00"
tasks = "24"
memory = "2GB"


# Submit many instances of main.py to different nodes
for ijk in range(n_sections):
    output = f"{routine}_output{ijk}.slurm"

    sbatch = open("par_submit.SBATCH","w")
    sbatch.write("#!/bin/bash \n")
    sbatch.write(f"#SBATCH {node} \n")
    sbatch.write(f"#SBATCH -J {ijk}_{routine} \n")
    sbatch.write(f"#SBATCH -o {output} \n")
    sbatch.write(f"#SBATCH -t {time} \n")
    # sbatch.write( "#SBATCH -N 1 \n")
    sbatch.write(f"#SBATCH --ntasks={tasks} \n")
    sbatch.write(f"#SBATCH --mem-per-cpu={memory} \n \n")

    # sbatch.write( "export OMP_NUM_THREADS=1 \n")
    # sbatch.write( "export MKL_NUM_THREADS=1 \n \n")

    # sbatch.write(f"cargo build --release \n \n")
    sbatch.write(f"/scratch/mtayl29/rad-in-rust/target/release/rad-in-rust {routine} {wc_norm} {log_g_bounds[ijk]} {log_g_bounds[ijk+1]} {ng} {nk} {n_kappa} {nf}")

    sbatch.close()

    print(f"g = {np.round(10 ** log_g_bounds[ijk], 3)} -> {np.round(10 ** log_g_bounds[ijk+1], 3)}")

    # Submit sbatch
    os.system("sbatch par_submit.SBATCH")