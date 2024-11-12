#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=48:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=256G                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=2                                        # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=20-25
#SBATCH --output=slurm_output/slurm-%A_%a-MafftAssembly.out

## Participants (0-indexed) 1 and 3 need more than 64G memory otherwise OOM crash, and 1 takes long than 1 hour. 128G/12hr worked.
## For 0,2 can use 64 GB and 2 cores
## For 1,3 apparently need 256G 12 hr and 1 core.
## For 4-12 can use multi cores and 32 GB memory.
## For 14,19 can use 32 GB and 2 cores.
## 

## Batch 1 = 0-2
## Podar = 3
## Tian = 4-10
## Ishiura = 11-12
## Batch 2 = 13-19
## Batch 3/4 = 20-25



## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load Python/3.11.3-gimkl-2022a
module load MAFFT/7.505-gimkl-2022a-with-extensions


## EXTRACT PARTICIPANT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="20JK4316D 21LE9817N 22NZ3355R podar_2023_21073LRa006 tian_2019_f1_iv_7 tian_2019_f1_iv_15 tian_2019_f2_ii_3 tian_2019_f4_ii_2 tian_2019_f5_ii_1 tian_2019_f5_ii_4 tian_2019_f9_ii_6 ishiura_2019_f9193_ii_5_raw ishiura_2019_f9193_ii_5_corrected GA5817C GQ1307H 20JK4322C 20JO4660C 22PP0924B 22RR5603Q 22RS5734T 311_merged 313_merged 317_merged 318_merged 319_merged 327_merged"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id

cd output/$participant_id
mkdir -p ../../consensuses

line_wrapping=75



## FLYE ASSEMBLY

if [ -s $participant_id"_5_long_reads_for_assembly.fastq" ]; then
    echo; echo \($(date +"%Y-%m-%d %T")\) assembling long allele via MAFFT alignment then custom consensus
    
    echo; echo \($(date +"%Y-%m-%d %T")\) long allele local alignment
    mafft --maxiterate 4000 --localpair $participant_id"_5_long_reads_for_assembly.fasta" > $participant_id"_6_long_mafft_alignment_local.msa"
    python ../../scripts/generate_custom_consensus.py -i $participant_id"_6_long_mafft_alignment_local.msa" -o $participant_id"_7_long_allele_consensus_local.fasta" -w $line_wrapping -v
    
    cp $participant_id"_7_long_allele_consensus_local.fasta" ../../consensuses
    
    echo; echo \($(date +"%Y-%m-%d %T")\) long allele global alignment
    mafft --maxiterate 4000 --globalpair $participant_id"_5_long_reads_for_assembly.fasta" > $participant_id"_6_long_mafft_alignment_global.msa"
    python ../../scripts/generate_custom_consensus.py -i $participant_id"_6_long_mafft_alignment_global.msa" -o $participant_id"_7_long_allele_consensus_global.fasta" -w $line_wrapping -v
    
    cp $participant_id"_7_long_allele_consensus_global.fasta" ../../consensuses
    
else
    echo; echo \($(date +"%Y-%m-%d %T")\) not assembling long allele \(no reads\)
fi


if [ -s $participant_id"_5_short_reads_for_assembly.fastq" ]; then
    echo; echo \($(date +"%Y-%m-%d %T")\) assembling short allele via MAFFT alignment then custom consensus
    
    echo; echo \($(date +"%Y-%m-%d %T")\) short allele local alignment
    mafft --maxiterate 4000 --localpair $participant_id"_5_short_reads_for_assembly.fasta" > $participant_id"_6_short_mafft_alignment_local.msa"
    python ../../scripts/generate_custom_consensus.py -i $participant_id"_6_short_mafft_alignment_local.msa" -o $participant_id"_7_short_allele_consensus_local.fasta" -w $line_wrapping -v
    
    cp $participant_id"_7_short_allele_consensus_local.fasta" ../../consensuses

    echo; echo \($(date +"%Y-%m-%d %T")\) short allele global alignment
    mafft --maxiterate 4000 --globalpair $participant_id"_5_short_reads_for_assembly.fasta" > $participant_id"_6_short_mafft_alignment_global.msa"
    python ../../scripts/generate_custom_consensus.py -i $participant_id"_6_short_mafft_alignment_global.msa" -o $participant_id"_7_short_allele_consensus_global.fasta" -w $line_wrapping -v
    
    cp $participant_id"_7_short_allele_consensus_global.fasta" ../../consensuses
    
else
    echo; echo \($(date +"%Y-%m-%d %T")\) not assembling short allele \(no reads\)
fi


if [ -s $participant_id"_5_hyper_reads_for_assembly.fastq" ]; then
    echo; echo \($(date +"%Y-%m-%d %T")\) assembling hyper allele via MAFFT alignment then custom consensus
    
    echo; echo \($(date +"%Y-%m-%d %T")\) hyper allele local alignment
    mafft --maxiterate 4000 --localpair $participant_id"_5_hyper_reads_for_assembly.fasta" > $participant_id"_6_hyper_mafft_alignment_local.msa"
    python ../../scripts/generate_custom_consensus.py -i $participant_id"_6_hyper_mafft_alignment_local.msa" -o $participant_id"_7_hyper_allele_consensus_local.fasta" -w $line_wrapping -v
    
    cp $participant_id"_7_hyper_allele_consensus_local.fasta" ../../consensuses

    echo; echo \($(date +"%Y-%m-%d %T")\) hyper allele global alignment
    mafft --maxiterate 4000 --globalpair $participant_id"_5_hyper_reads_for_assembly.fasta" > $participant_id"_6_hyper_mafft_alignment_global.msa"
    python ../../scripts/generate_custom_consensus.py -i $participant_id"_6_hyper_mafft_alignment_global.msa" -o $participant_id"_7_hyper_allele_consensus_global.fasta" -w $line_wrapping -v
    
    cp $participant_id"_7_hyper_allele_consensus_global.fasta" ../../consensuses
    
else
    echo; echo \($(date +"%Y-%m-%d %T")\) not assembling hyper allele \(no reads\)
fi


echo; echo \($(date +"%Y-%m-%d %T")\) done
