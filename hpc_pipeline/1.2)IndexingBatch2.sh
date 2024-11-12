#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=48:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=48000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-6
#SBATCH --output=slurm_output/slurm-%A_%a-IndexingBatch2.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load minimap2/2.27-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0


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




cd output

mkdir -p $participant_id
cd $participant_id


echo \($(date +"%Y-%m-%d %T")\) Mapping participant reads to reference genome
minimap2 -a -y -t $OMP_NUM_THREADS -x map-ont ../hg38.mmi ../../input/emma/batch_2/$experiment_id"_pass_sup.fastq.gz" > $participant_id"_1_all_reads.sam"

echo \($(date +"%Y-%m-%d %T")\) Sorting and binarising aligned reads
samtools sort -@ $OMP_NUM_THREADS $participant_id"_1_all_reads.sam" > $participant_id"_1_all_reads.bam"

echo \($(date +"%Y-%m-%d %T")\) Indexing sorted reads
samtools index -@ $OMP_NUM_THREADS $participant_id"_1_all_reads.bam"

echo \($(date +"%Y-%m-%d %T")\) Deleting uncompressed reads
rm $participant_id"_1_all_reads.sam"

echo; echo \($(date +"%Y-%m-%d %T")\) done
