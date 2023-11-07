#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --output=output/array_%A_%a.out
#SBATCH --array=1-50

export FILE="data_batch_4_20_0.1" # "data_batch_4_20_0.3"
# "data_batch_6_10_0.1" # "data_batch_6_10_0.3"
# "data_batch_6_15_0.1" # "data_batch_6_15_0.3"
# "data_batch_6_20_0.1" # "data_batch_6_20_0.3"
# "data_batch_8_10_0.1" # "data_batch_8_10_0.3"
# "data_batch_8_15_0.1" # "data_batch_8_15_0.3"
# "data_batch_8_20_0.1" # "data_batch_8_20_0.3"
export K=1  # k 1 2 3
export L=$SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
mkdir results/$FILE

module load julia
srun julia --project=. arraymain.jl $K $L $FILE
awk '(NR == 1) || (FNR > 1)' results/$FILE/results_${FILE}_k${K}_*.csv > results/$FILE/combined_results_${FILE}_k$K.csv