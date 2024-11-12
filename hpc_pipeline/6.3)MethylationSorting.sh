#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=24:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-15
#SBATCH --output=slurm_output/slurm-%A_%a-MethylationSorting.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load SAMtools/1.19-GCC-12.3.0
module load minimap2/2.27-GCC-12.3.0


## SET PARTICIPANT ID
i=$SLURM_ARRAY_TASK_ID
participants="20JK4316D 21LE9817N 22NZ3355R GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T 311 313 317 318 319 327"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id; echo

cd dorado/$participant_id



echo; echo \($(date +"%Y-%m-%d %T")\) generating list of all read ids from initial basecalling
samtools fastq -@ $OMP_NUM_THREADS -T MM,ML $participant_id"_2_calls.bam" > $participant_id"_2_calls.fastq" 
awk "(NR-1) % 4 == 0" $participant_id"_2_calls.fastq" > $participant_id"_2_calls_read_ids.txt"


echo; echo \($(date +"%Y-%m-%d %T")\) mapping and sorting reads
minimap2 -a -y -x map-ont ../../output/hg38.mmi $participant_id"_2_calls.fastq" > $participant_id"_3_calls_indexed.sam"
samtools sort -@ $OMP_NUM_THREADS $participant_id"_3_calls_indexed.sam" > $participant_id"_3_calls_sorted.bam"
samtools index -@ $OMP_NUM_THREADS $participant_id"_3_calls_sorted.bam"
rm $participant_id"_3_calls_indexed.sam"


echo \($(date +"%Y-%m-%d %T")\) done
