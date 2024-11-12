#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=12:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)
#SBATCH --partition=large
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-2
#SBATCH --output=slurm_output/slurm-%A_%a-IndexingBatch1.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load minimap2/2.27-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0


## EXTRACT PARTICIPANT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="20JK4316D 21LE9817N 22NZ3355R"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id; echo





cd output

mkdir -p $participant_id
cd $participant_id



echo \($(date +"%Y-%m-%d %T")\) Mapping participant reads to reference genome
minimap2 -a -x map-ont ../hg38.mmi ../../input/emma/batch_1/$participant_id"_pass.fastq.gz" > $participant_id"_1_all_reads.sam"

echo \($(date +"%Y-%m-%d %T")\) Sorting and binarising aligned reads
samtools sort $participant_id"_1_all_reads.sam" > $participant_id"_1_all_reads.bam"

echo \($(date +"%Y-%m-%d %T")\) Indexing sorted reads
samtools index $participant_id"_1_all_reads.bam"

echo; echo \($(date +"%Y-%m-%d %T")\) done