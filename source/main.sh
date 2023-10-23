#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --output=output.out

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
module purge
module load julia
julia --project=. main.jl