#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=48:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=32000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)
#SBATCH --gpus-per-node=A100:1

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-13
#SBATCH --output=slurm_output/slurm-%A_%a-Dorado.out



## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load Dorado/0.7.3
module load CUDA/12.2.2


## SET EXPERIMENT ID
i=$SLURM_ARRAY_TASK_ID
experiments="PGAXHC230412 PGAXHC240012 PGAXHX240013 GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T PGAXOW240379 PGAXOW240380 PGAXOW240412 PGAXOW240413"
experiments=($experiments)
experiment_id=${experiments[$i]}
echo Experiment ID: $experiment_id; echo


cd dorado/$experiment_id

dorado basecaller sup,5mCG_5hmCG $experiment_id"_1_reads.pod5" --kit-name SQK-NBD114-24 > $experiment_id"_2_calls.bam"


module load SAMtools/1.19-GCC-12.3.0
samtools fastq -@ $OMP_NUM_THREADS -T MM,ML $experiment_id"_2_calls.bam" > $experiment_id"_2_calls.fastq"
awk "(NR-1) % 4 == 0" $experiment_id"_2_calls.fastq" > $experiment_id"_2_calls_read_ids.txt"

echo \($(date +"%Y-%m-%d %T")\) done
