#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --output=array_%A_%a.out
#SBATCH --array=0-5

export FILE = "test_batch.txt"
export K = 2
export L = $SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

module load julia
srun julia --project=. arraymain.jl $K $L $FILE
awk '(NR == 1) || (FNR > 1)' results/$FILE/*.csv > combined.csv