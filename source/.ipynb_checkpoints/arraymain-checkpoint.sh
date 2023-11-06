#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=8
#SBATCH --output=output/array_%A_%a.out
#SBATCH --array=1-5

export FILE="test_batch"
export K=1
export L=$SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
mkdir results/$FILE

module load julia
srun julia --project=. arraymain.jl $K $L $FILE
awk '(NR == 1) || (FNR > 1)' results/$FILE/results_${FILE}_k${K}_*.csv > results/$FILE/combined_results_${FILE}_k$K.csv