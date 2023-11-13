#!/bin/bash
#SBATCH --time=00:1:00
#SBATCH --mem=1G

# export FILE="data_batch_4_10_0.1" 
# export FILE="data_batch_4_10_0.3"
# export FILE="data_batch_4_15_0.1" 
# export FILE="data_batch_4_15_0.3"
# export FILE="data_batch_4_20_0.1" 
# export FILE="data_batch_4_20_0.3"
# n = 6
# export FILE="data_batch_6_10_0.1" 
# export FILE="data_batch_6_10_0.3"
# export FILE="data_batch_6_15_0.1" 
# export FILE="data_batch_6_15_0.3"
# export FILE="data_batch_6_20_0.1" 
# export FILE="data_batch_6_20_0.3"
# n = 8
# export FILE="data_batch_8_10_0.1" 
# export FILE="data_batch_8_10_0.3"
# export FILE="data_batch_8_15_0.1" 
# export FILE="data_batch_8_15_0.3"
# export FILE="data_batch_8_20_0.1" 
# 
export FILE="data_batch_8_20_0.3"
# k
export K=2
awk '(NR == 1) || (FNR > 1)' results/$FILE/individual/results_${FILE}_k${K}_*.csv > results/$FILE/combined_results_${FILE}_k$K.csv
