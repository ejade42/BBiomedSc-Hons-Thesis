#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=00:30:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-5
#SBATCH --output=slurm_output/slurm-%A_%a-BatchMerging.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a
module load SAMtools/1.19-GCC-12.3.0
module load seqtk/1.4-GCC-11.3.0



## EXTRACT PARTICIPANT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="311 313 317 318 319 327"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id


cd output
mkdir -p $participant_id"_merged"

files=("_3_NOTCH2NLC_primary_reads.fastq"
       "_4_NOTCH2NLC_clipped_reads_forward.fastq"
       "_4_NOTCH2NLC_clipped_reads_reverse.fastq")

for filename in ${files[*]}; 
do 
    echo Merging for file: $filename
    cat $participant_id"/"$participant_id$filename $participant_id"_2/"$participant_id"_2"$filename > $participant_id"_merged/"$participant_id"_merged"$filename
done

echo \($(date +"%Y-%m-%d %T")\) done
