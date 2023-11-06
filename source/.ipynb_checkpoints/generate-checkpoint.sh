#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --output=output.out

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
module load julia
julia --project=. generate_ins.jl