#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=00:30:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-5
#SBATCH --output=slurm_output/slurm-%A_%a-PostMergeFiltering.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a
module load SAMtools/1.19-GCC-12.3.0
module load seqtk/1.4-GCC-11.3.0



## EXTRACT PARTICIPANT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="311_merged 313_merged 317_merged 318_merged 319_merged 327_merged"
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