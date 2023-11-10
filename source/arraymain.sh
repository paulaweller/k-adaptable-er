#!/bin/bash
#SBATCH --time=02:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --output=output/array_%A_%a.out
#SBATCH --array=1-50
#SBATCH --exclude=milan14

# export FILE="data_batch_4_10_0.1" 
# export FILE="data_batch_4_10_0.3"
# export FILE="data_batch_4_15_0.1" 
# export FILE="data_batch_4_15_0.3"
# export FILE="data_batch_4_20_0.1" 
# export FILE="data_batch_4_20_0.3"
# n = 6
# export FILE="data_batch_6_10_0.1" 
# export FILE="data_batch_6_10_0.3"
#  export FILE="data_batch_6_15_0.1" 
# export FILE="data_batch_6_15_0.3"
# export FILE="data_batch_6_20_0.1" 
# export FILE="data_batch_6_20_0.3"
# n = 8
# export FILE="data_batch_8_10_0.1" 
# export FILE="data_batch_8_10_0.3"
# export FILE="data_batch_8_15_0.1" 
# export FILE="data_batch_8_15_0.3"
# export FILE="data_batch_8_20_0.1" 
# export FILE="data_batch_8_20_0.3"
# k
export K=3
export L=$SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
mkdir results/$FILE
mkdir results/$FILE/individual

module load julia
srun julia --project=. arraymain.jl $K $L $FILE
awk '(NR == 1) || (FNR > 1)' results/$FILE/individual/results_${FILE}_k${K}_*.csv > results/$FILE/combined_results_${FILE}_k$K.csv

# rm results/data_batch_4_10_0.1/individual/results_data_batch_4_10_0.1_k1_*.csv
####### find . -type f -name results_data_batch_6_10_0.1_\* -exec rm {} \;