#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=24:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=8000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=2                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-6
#SBATCH --output=slurm_output/slurm-%A_%a-Pod5ProcessingBatch2.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
source ./blue-crab-venv/bin/activate


## EXTRACT PARTICIPANT AND EXPERIMENT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id

experiments="PGXXOW240233 PGXXOW240234 PGXXOW240235 PGXXOW240236 PGXXOW240237 PGXXOW240238 PGXXOW240239"
experiments=($experiments)
experiment_id=${experiments[$i]}
echo Experiment ID: $experiment_id; echo


cd dorado
mkdir -p $participant_id
cd $participant_id


echo \($(date +"%Y-%m-%d %T")\) Using blue-crab to convert blow5 to pod5
rm -f $participant_id"_1_reads.pod5"
blue-crab s2p ../../input/emma/batch_2/$experiment_id"_reads.blow5" -o $participant_id"_1_reads.pod5"


echo \($(date +"%Y-%m-%d %T")\) done
