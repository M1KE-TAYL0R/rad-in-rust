import os
import numpy as np

# Simulation parameters
n_sections = 4
g_min_log =  1
g_max_log =  2
log_g_bounds = np.linspace(g_min_log,g_max_log,n_sections+1)
ng = 4
nf = 2
nk = 512
n_kappa = 501

print(f"Running main.py with {ng*n_sections} couplings and {nk} k-points")

# sbatch parameters
node = "-p action -A action"
# node = "-p standard"
job = "rad_rust"
time = "2-00:00:00"
tasks = "24"
memory = "60GB"


# Submit many instances of main.py to different nodes
for ijk in range(n_sections):
    output = f"output{ijk}.slurm"

    sbatch = open("SUBMIT_PAR.SBATCH","w")
    sbatch.write("#!/bin/bash \n")
    sbatch.write(f"#SBATCH {node} \n")
    sbatch.write(f"#SBATCH -J {ijk}_{job} \n")
    sbatch.write(f"#SBATCH -o {output} \n")
    sbatch.write(f"#SBATCH -t {time} \n")
    sbatch.write( "#SBATCH -N 1 \n")
    sbatch.write(f"#SBATCH --ntasks-per-node={tasks} \n")
    sbatch.write(f"#SBATCH --mem {memory} \n \n")

    # sbatch.write( "export OMP_NUM_THREADS=1 \n")
    # sbatch.write( "export MKL_NUM_THREADS=1 \n \n")

    sbatch.write("cargo build --release")
    sbatch.write(f"./target/release/rad_in_rust {log_g_bounds[ijk]} {log_g_bounds[ijk+1]} {ng} {nk} {n_kappa} {nf}")

    sbatch.close()

    print(f"g = {np.round(10 ** log_g_bounds[ijk], 3)} -> {np.round(10 ** log_g_bounds[ijk+1], 3)}")

    # Submit sbatch
    os.system("sbatch par_submit.SBATCH")