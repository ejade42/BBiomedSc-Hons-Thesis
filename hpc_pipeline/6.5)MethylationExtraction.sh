#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=1:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=8000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=2                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-15
#SBATCH --output=slurm_output/slurm-%A_%a-MethylationExtraction.out



## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load modkit/0.2.5-GCC-12.3.0
module load Python/3.11.6-foss-2023a
module load SAMtools/1.19-GCC-12.3.0


## SET PARTICIPANT ID
i=$SLURM_ARRAY_TASK_ID
participants="20JK4316D 21LE9817N 22NZ3355R GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T 311 313 317 318 319 327"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id; echo

cd dorado/$participant_id


## Process full-length reads
echo \($(date +"%Y-%m-%d %T")\) extracting methylation information for full-length reads
python3 ../../scripts/parse_modifications_fastq.py -i $participant_id"_4_NOTCH2NLC_reads.fastq" -o-r $participant_id"_4_NOTCH2NLC_reads_information.csv" -o-m $participant_id"_4_NOTCH2NLC_reads_methylation.txt" -v




## Process repeat-region trimmed reads
echo; echo \($(date +"%Y-%m-%d %T")\) extracting methylation information for repeat region trimmed reads

samtools view -F 16 --threads $OMP_NUM_THREADS $participant_id"_6_NOTCH2NLC_trimmed_reads.bam" -u | \
samtools fastq - | \
awk "(NR-1) % 4 == 0" - > $participant_id"_6_NOTCH2NLC_trimmed_reads_forward_ids.txt"

samtools view -f 16 --threads $OMP_NUM_THREADS $participant_id"_6_NOTCH2NLC_trimmed_reads.bam" -u | \
samtools fastq - | \
awk "(NR-1) % 4 == 0" - > $participant_id"_6_NOTCH2NLC_trimmed_reads_reverse_ids.txt"


python3 ../../scripts/parse_modifications_fastq.py -i $participant_id"_6_NOTCH2NLC_trimmed_reads.fastq" -o-r $participant_id"_6_NOTCH2NLC_trimmed_reads_information.csv" -o-m $participant_id"_6_NOTCH2NLC_trimmed_reads_methylation.txt" -v


module load R/4.2.1-gimkl-2022a
if [ -s $participant_id"_6_NOTCH2NLC_trimmed_reads_methylation.txt" ]; then
    echo; echo \($(date +"%Y-%m-%d %T")\) exporting methylation to csv
    Rscript ../../scripts/extract_methylation.R $participant_id $participant_id"_6_NOTCH2NLC_trimmed_reads" $participant_id"_methylation_information.csv"

else
    echo; echo \($(date +"%Y-%m-%d %T")\) not exporting methylation information \(no reads\)
    echo "participant_id","read_id","direction","length","sequence","quality","modification_types","C_hm_locations","C_hm_probabilities","C_m_locations","C_m_probabilities" > $participant_id"_methylation_information.csv"
fi

echo; echo \($(date +"%Y-%m-%d %T")\) done
