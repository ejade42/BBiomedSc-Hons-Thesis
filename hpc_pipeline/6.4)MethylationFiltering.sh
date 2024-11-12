#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=8:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=12000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-15
#SBATCH --output=slurm_output/slurm-%A_%a-MethylationFiltering.out



## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load SAMtools/1.19-GCC-12.3.0
module load modkit/0.2.5-GCC-12.3.0


## SET PARTICIPANT ID
i=$SLURM_ARRAY_TASK_ID
participants="20JK4316D 21LE9817N 22NZ3355R GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T 311 313 317 318 319 327"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id; echo

cd dorado/$participant_id



## SUBSET TO NOTCH2NLC PRIMARY MAPPINGS
echo; echo \($(date +"%Y-%m-%d %T")\) subsetting reads
samtools view -L ../../input/NOTCH2NLC_promoters.bed -F 256 -q 10 --threads $OMP_NUM_THREADS $participant_id"_3_calls_sorted.bam" -o $participant_id"_4_NOTCH2NLC_reads.bam"
samtools index -@ $OMP_NUM_THREADS $participant_id"_4_NOTCH2NLC_reads.bam"
samtools fastq -@ $OMP_NUM_THREADS $participant_id"_4_NOTCH2NLC_reads.bam" -T MM,ML > $participant_id"_4_NOTCH2NLC_reads.fastq"



## SUBSET FURTHER TO ONLY READS COMPLETELY OVERLAPPING THE REPEAT, AND TRIM READS TO JUST THE REPEAT
echo; echo \($(date +"%Y-%m-%d %T")\) trimming reads
samtools view -@ $OMP_NUM_THREADS -L ../../input/NOTCH2NLC_region_start.bed $participant_id"_4_NOTCH2NLC_reads.bam" -u | \
samtools view -@ $OMP_NUM_THREADS -L ../../input/NOTCH2NLC_region_end.bed -u - |
samtools ampliconclip --both-ends --hard-clip -b ../../input/NOTCH2NLC_region_to_clip.bed -o $participant_id"_5_1_NOTCH2NLC_clipped_reads.bam" -

## REPAIR TAGS
echo; echo \($(date +"%Y-%m-%d %T")\) repairing tags
samtools sort -@ $OMP_NUM_THREADS -n $participant_id"_4_NOTCH2NLC_reads.bam" > $participant_id"_5_2_name_sorted_NOTCH2NLC_reads.bam"
samtools sort -@ $OMP_NUM_THREADS -n $participant_id"_5_1_NOTCH2NLC_clipped_reads.bam" > $participant_id"_5_2_name_sorted_clipped_reads.bam"

modkit repair -t $OMP_NUM_THREADS -d $participant_id"_5_2_name_sorted_NOTCH2NLC_reads.bam" -a $participant_id"_5_2_name_sorted_clipped_reads.bam" -o $participant_id"_5_3_name_sorted_repaired_reads.bam" --log-filepath "modkit_log.txt"

samtools sort -@ $OMP_NUM_THREADS $participant_id"_5_3_name_sorted_repaired_reads.bam" > $participant_id"_6_NOTCH2NLC_trimmed_reads.bam"
samtools fastq -@ $OMP_NUM_THREADS $participant_id"_6_NOTCH2NLC_trimmed_reads.bam" -T MM,ML > $participant_id"_6_NOTCH2NLC_trimmed_reads.fastq"


echo; echo \($(date +"%Y-%m-%d %T")\) done

