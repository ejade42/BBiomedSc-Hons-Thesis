#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=00:05:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=4000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=4                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END

#SBATCH --output=slurm_output/slurm-%j-ConsensusProcessing.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a


## Batch 1 participants
batch_1_participants="20JK4316D 21LE9817N 22NZ3355R podar_2023_21073LRa006 tian_2019_f1_iv_7 tian_2019_f1_iv_15 tian_2019_f2_ii_3 tian_2019_f4_ii_2 tian_2019_f5_ii_1 tian_2019_f5_ii_4 tian_2019_f9_ii_6 ishiura_2019_f9193_ii_5_raw ishiura_2019_f9193_ii_5_corrected"

## Batch 2 participants
batch_2_participants="GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T"

## Batch 3 participants
batch_3_participants="311 313 317 318 319 327"

## Batch 4 participants
batch_4_participants="311_2 313_2 317_2 318_2 319_2 327_2"

## Batch 3-4 merged participants
batch_3_4_participants="311_merged 313_merged 317_merged 318_merged 319_merged 327_merged"

cd output


## Run read information collation script
echo; echo \($(date +"%Y-%m-%d %T")\) collating read information
Rscript ../scripts/collate_read_information.R $batch_1_participants "batch_1_read_information_all_participants.csv"
Rscript ../scripts/collate_read_information.R $batch_2_participants "batch_2_read_information_all_participants.csv"
Rscript ../scripts/collate_read_information.R $batch_3_participants "batch_3_read_information_all_participants.csv"
Rscript ../scripts/collate_read_information.R $batch_4_participants "batch_4_read_information_all_participants.csv"
Rscript ../scripts/collate_read_information.R $batch_3_4_participants "batch_3_4_read_information_all_participants.csv"

cd ..


echo; echo \($(date +"%Y-%m-%d %T")\) done