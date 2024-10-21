#!/bin/bash
#SBATCH --time=04:20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=output_rio/rio_array_%A_%a.out
#SBATCH --array=1-5

# supply product: Food, Water, Hygiene, Cleaning, Mattress, Medicine
export PRODUCT="Medicine"
# k
export K=$SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
mkdir results/rio/$PRODUCT
mkdir results/rio/$PRODUCT/individual

module load julia
srun julia --project=. rio_main.jl $K $PRODUCT

# rm results/data_batch_4_10_0.1/individual/results_data_batch_4_10_0.1_k1_*.csv
####### find . -type f -name results_data_batch_6_10_0.1_\* -exec rm {} \;