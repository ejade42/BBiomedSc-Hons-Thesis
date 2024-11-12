#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=1:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=4                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-1
#SBATCH --output=slurm_output/slurm-%A_%a-IndexingIshiura.out



## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load minimap2/2.27-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0


## EXTRACT PARTICIPANT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="ishiura_2019_f9193_ii_5_raw ishiura_2019_f9193_ii_5_corrected"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id; echo



cd output

mkdir -p $participant_id
cd $participant_id


echo \($(date +"%Y-%m-%d %T")\) Mapping participant reads to reference genome - change "raw" or "appended" to select input
minimap2 -a -t $OMP_NUM_THREADS -x map-ont ../hg38.mmi ../../input/ishiura_2019/appended/$participant_id".fasta" > $participant_id"_1_all_reads.sam"


echo \($(date +"%Y-%m-%d %T")\) Sorting and binarising aligned reads
samtools sort -@ $OMP_NUM_THREADS $participant_id"_1_all_reads.sam" > $participant_id"_1_all_reads.bam"


echo \($(date +"%Y-%m-%d %T")\) Indexing sorted reads
samtools index $participant_id"_1_all_reads.bam"


echo; echo \($(date +"%Y-%m-%d %T")\) done
