#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=6:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=8000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-6
#SBATCH --output=slurm_output/slurm-%A_%a-Demuxing.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load SAMtools/1.19-GCC-12.3.0
module load Dorado/0.7.3
module load CUDA/12.2.2


i=$SLURM_ARRAY_TASK_ID
experiments="PGAXHC230412 PGAXHC240012 PGAXHX240013 PGAXOW240379 PGAXOW240380 PGAXOW240412 PGAXOW240413"
experiments=($experiments)
experiment_id=${experiments[$i]}
echo Experiment ID: $experiment_id; echo

cd dorado/$experiment_id

dorado demux --output-dir $experiment_id"_demuxed/" --no-classify $experiment_id"_2_calls.bam"


echo \($(date +"%Y-%m-%d %T")\) done
