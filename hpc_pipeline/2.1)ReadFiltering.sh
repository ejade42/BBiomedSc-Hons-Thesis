#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=00:30:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-31
#SBATCH --output=slurm_output/slurm-%A_%a-Filtering.out

## Batch 1 = 0-2
## Podar = 3
## Tian = 4-10
## Ishiura = 11-12
## Batch 2 = 13-19
## Batch 3 = 20-25


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a
module load SAMtools/1.19-GCC-12.3.0
module load seqtk/1.4-GCC-11.3.0



## EXTRACT PARTICIPANT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="20JK4316D 21LE9817N 22NZ3355R podar_2023_21073LRa006 tian_2019_f1_iv_7 tian_2019_f1_iv_15 tian_2019_f2_ii_3 tian_2019_f4_ii_2 tian_2019_f5_ii_1 tian_2019_f5_ii_4 tian_2019_f9_ii_6 ishiura_2019_f9193_ii_5_raw ishiura_2019_f9193_ii_5_corrected GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T 311 313 317 318 319 327 311_2 313_2 317_2 318_2 319_2 327_2"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id


cd output/$participant_id


## STEP 1: Filter to only reads overlapping both ends of the NOTCH2NLC repeat
## The A of the ATG is at chr1:149,390,788
## The first GGC in HG38 is at chr1:149,390,803-149,390,805
## The last GGC in HG38 is at chr1:149,390,839-149,390,841
## If we take a little bit more downstream, we could cut at the first T at chr1:149,390,853
## So, let's try to get everything overlapping with 788 AND with 853.


echo; echo \($(date +"%Y-%m-%d %T")\) filtering to reads completely overlapping the repeat only
samtools view -L ../../input/NOTCH2NLC_region_start.bed --threads $OMP_NUM_THREADS $participant_id"_1_all_reads.bam" -u | \
samtools view -L ../../input/NOTCH2NLC_region_end.bed --threads $OMP_NUM_THREADS -o $participant_id"_2_NOTCH2NLC_reads.bam" -






## STEP 2: Filter to only primary reads (avoid NOTCH2NLA/B)
## 256 means "not primary", so use -F to exclude all non-primaries

echo; echo \($(date +"%Y-%m-%d %T")\) filtering to primary reads only
samtools view -F 256 --threads $OMP_NUM_THREADS -o $participant_id"_3_NOTCH2NLC_primary_reads.bam" $participant_id"_2_NOTCH2NLC_reads.bam"
samtools index $participant_id"_3_NOTCH2NLC_primary_reads.bam"
samtools fastq $participant_id"_3_NOTCH2NLC_primary_reads.bam" > $participant_id"_3_NOTCH2NLC_primary_reads.fastq"





## STEP 3: Clip reads to only the selected NOTCH2NLC repeat region

echo; echo \($(date +"%Y-%m-%d %T")\) clipping reads
samtools ampliconclip --both-ends --hard-clip -b ../../input/NOTCH2NLC_region_to_clip.bed $participant_id"_3_NOTCH2NLC_primary_reads.bam" > $participant_id"_4_NOTCH2NLC_clipped_reads.bam"
samtools index $participant_id"_4_NOTCH2NLC_clipped_reads.bam"
samtools fastq $participant_id"_4_NOTCH2NLC_clipped_reads.bam" > $participant_id"_4_NOTCH2NLC_clipped_reads.fastq"






## STEP 4: Split clipped reads into forward and reverse

samtools view -F 16 --threads $OMP_NUM_THREADS $participant_id"_4_NOTCH2NLC_clipped_reads.bam" -u | \
samtools fastq - > $participant_id"_4_NOTCH2NLC_clipped_reads_forward.fastq"

samtools view -f 16 --threads $OMP_NUM_THREADS $participant_id"_4_NOTCH2NLC_clipped_reads.bam" -u | \
samtools fastq - > $participant_id"_4_NOTCH2NLC_clipped_reads_reverse.fastq"






## STEP 5: Reformat reads (via R)

echo; echo \($(date +"%Y-%m-%d %T")\) reformatting into csv
Rscript ../../scripts/make_csv_from_fastq.R $participant_id 200 800
Rscript ../../scripts/make_list_of_long_short_reads.R $participant_id



## STEP 6: Extract reads for assembly

echo; echo \($(date +"%Y-%m-%d %T")\) extracting long and short reads for assembly
seqtk subseq $participant_id"_3_NOTCH2NLC_primary_reads.fastq" $participant_id"_long_reads.txt" > $participant_id"_5_long_reads_for_assembly.fastq"
seqtk subseq $participant_id"_3_NOTCH2NLC_primary_reads.fastq" $participant_id"_short_reads.txt" > $participant_id"_5_short_reads_for_assembly.fastq"
seqtk subseq $participant_id"_3_NOTCH2NLC_primary_reads.fastq" $participant_id"_hyper_reads.txt" > $participant_id"_5_hyper_reads_for_assembly.fastq"


seqtk seq -a $participant_id"_5_long_reads_for_assembly.fastq" > $participant_id"_5_long_reads_for_assembly.fasta"
seqtk seq -a $participant_id"_5_short_reads_for_assembly.fastq" > $participant_id"_5_short_reads_for_assembly.fasta"
seqtk seq -a $participant_id"_5_hyper_reads_for_assembly.fastq" > $participant_id"_5_hyper_reads_for_assembly.fasta"




echo; echo \($(date +"%Y-%m-%d %T")\) done