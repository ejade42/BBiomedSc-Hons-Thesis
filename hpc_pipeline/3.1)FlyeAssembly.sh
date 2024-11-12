#!/bin/bash
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=01:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=32000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-25
#SBATCH --output=slurm_output/slurm-%A_%a-FlyeAssembly.out

## Batch 1 = 0-2
## Podar = 3
## Tian = 4-10
## Ishiura = 11-12
## Batch 2 = 13-19
## Batch 3/4 = 20-25



## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a
module load bzip2/1.0.8-GCCcore-11.3.0
module load libxml2/2.9.10-GCCcore-11.3.0
module load Flye/2.9.3-gimkl-2022a-Python-3.11.3


## EXTRACT PARTICIPANT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="20JK4316D 21LE9817N 22NZ3355R podar_2023_21073LRa006 tian_2019_f1_iv_7 tian_2019_f1_iv_15 tian_2019_f2_ii_3 tian_2019_f4_ii_2 tian_2019_f5_ii_1 tian_2019_f5_ii_4 tian_2019_f9_ii_6 ishiura_2019_f9193_ii_5_raw ishiura_2019_f9193_ii_5_corrected GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T 311_merged 313_merged 317_merged 318_merged 319_merged 327_merged"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id

cd output/$participant_id
mkdir -p ../../consensuses



## FLYE ASSEMBLY

if [ -s $participant_id"_5_long_reads_for_assembly.fastq" ]; then
    echo; echo \($(date +"%Y-%m-%d %T")\) assembling long allele via Flye
    flye --nano-hq $participant_id"_5_long_reads_for_assembly.fastq" --out-dir $participant_id"_assembly_long_allele" --threads $OMP_NUM_THREADS --min-overlap 1000
    
    cp $participant_id"_assembly_long_allele/assembly.fasta" $participant_id"_7_long_allele_flye.fasta"
    cp $participant_id"_7_long_allele_flye.fasta" ../../consensuses

else
    echo; echo \($(date +"%Y-%m-%d %T")\) not assembling long allele via Flye as no reads provided
fi


if [ -s $participant_id"_5_short_reads_for_assembly.fastq" ]; then
    echo; echo \($(date +"%Y-%m-%d %T")\) assembling short allele via Flye
    flye --nano-hq $participant_id"_5_short_reads_for_assembly.fastq" --out-dir $participant_id"_assembly_short_allele" --threads $OMP_NUM_THREADS --min-overlap 1000
    
    cp $participant_id"_assembly_short_allele/assembly.fasta" $participant_id"_7_short_allele_flye.fasta"
    cp $participant_id"_7_short_allele_flye.fasta" ../../consensuses

else
    echo; echo \($(date +"%Y-%m-%d %T")\) not assembling short allele via Flye as no reads provided
fi


if [ -s $participant_id"_5_hyper_reads_for_assembly.fastq" ]; then
    echo; echo \($(date +"%Y-%m-%d %T")\) assembling hyper allele via Flye
    flye --nano-hq $participant_id"_5_hyper_reads_for_assembly.fastq" --out-dir $participant_id"_assembly_hyper_allele" --threads $OMP_NUM_THREADS --min-overlap 1000
    
    cp $participant_id"_assembly_hyper_allele/assembly.fasta" $participant_id"_7_hyper_allele_flye.fasta"
    cp $participant_id"_7_hyper_allele_flye.fasta" ../../consensuses

else
    echo; echo \($(date +"%Y-%m-%d %T")\) not assembling hyper allele via Flye as no reads provided
fi


echo; echo \($(date +"%Y-%m-%d %T")\) done
