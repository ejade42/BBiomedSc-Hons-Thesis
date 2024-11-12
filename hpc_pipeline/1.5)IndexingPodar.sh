#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=48:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=40000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END

#SBATCH --output=slurm_output/slurm-%j-IndexingPodar.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load minimap2/2.27-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0


## EXTRACT PARTICIPANT TO USE
participant_id=podar_2023_21073LRa006
echo Participant ID: $participant_id; echo





cd output

mkdir -p $participant_id
cd $participant_id



echo \($(date +"%Y-%m-%d %T")\) Copying provided bam to output directory as fastq
samtools fastq ../../input/podar_2023/21073LRa006_01.bam > $participant_id"_1_all_reads.fastq"


echo \($(date +"%Y-%m-%d %T")\) Mapping participant reads to reference genome
minimap2 -a -t $OMP_NUM_THREADS -x map-ont ../hg38.mmi $participant_id"_1_all_reads.fastq" > $participant_id"_1_all_reads.sam"


echo \($(date +"%Y-%m-%d %T")\) Sorting and binarising aligned reads
samtools sort -@ $OMP_NUM_THREADS $participant_id"_1_all_reads.sam" > $participant_id"_1_all_reads.bam"


echo \($(date +"%Y-%m-%d %T")\) Indexing sorted reads
samtools index $participant_id"_1_all_reads.bam"


echo; echo \($(date +"%Y-%m-%d %T")\) done