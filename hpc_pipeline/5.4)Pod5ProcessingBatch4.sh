#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=6:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=8000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=2                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-1
#SBATCH --output=slurm_output/slurm-%A_%a-Pod5ProcessingBatch4.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
source ./blue-crab-venv/bin/activate


## SET EXPERIMENT ID
i=$SLURM_ARRAY_TASK_ID
experiments="PGAXOW240412 PGAXOW240413"
experiments=($experiments)
experiment_id=${experiments[$i]}
echo Experiment ID: $experiment_id; echo

cd dorado
mkdir -p $experiment_id
cd $experiment_id


echo \($(date +"%Y-%m-%d %T")\) Using blue-crab to convert blow5 to pod5
rm -f $experiment_id"_reads.pod5"
blue-crab s2p ../../input/emma/batch_4/$experiment_id"_reads.blow5" -o $experiment_id"_1_reads.pod5"


echo \($(date +"%Y-%m-%d %T")\) done
