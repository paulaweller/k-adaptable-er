#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --output=outputnew.out
#SBATCH --cpus-per-task=8

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
module load julia
srun julia --project=. main.jl