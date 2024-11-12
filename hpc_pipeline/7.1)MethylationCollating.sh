#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=00:05:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=4000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=4                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END

#SBATCH --output=slurm_output/slurm-%j-MethylationCollating.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a


## Batch 1 participants
participants="20JK4316D 21LE9817N 22NZ3355R 311 313 317 318 319 327"

cd dorado

## Run read information collation script
echo; echo \($(date +"%Y-%m-%d %T")\) collating read information
Rscript ../scripts/collate_methylation_information.R $participants "methylation_information_all_participants.csv"

cd ..


echo; echo \($(date +"%Y-%m-%d %T")\) done