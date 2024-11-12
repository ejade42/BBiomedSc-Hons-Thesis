#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=24:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=24000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-5
#SBATCH --output=slurm_output/slurm-%A_%a-IndexingBatch3.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load minimap2/2.27-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0


## EXTRACT PARTICIPANT AND EXPERIMENT TO USE
i=$SLURM_ARRAY_TASK_ID
participants="311_2 313_2 317_2 318_2 319_2 327_2"
participants=($participants)
participant_id=${participants[$i]}
echo Participant ID: $participant_id

experiments="PGAXOW240412 PGAXOW240413 PGAXOW240413 PGAXOW240412 PGAXOW240413 PGAXOW240412"
experiments=($experiments)
experiment_id=${experiments[$i]}
echo Experiment ID: $experiment_id

barcodes="14 11 12 10 13 09"
barcodes=($barcodes)
barcode=${barcodes[$i]}
echo Barcode: $barcode




cd output

mkdir -p $participant_id
cd $participant_id


echo \($(date +"%Y-%m-%d %T")\) Mapping participant reads to reference genome
minimap2 -a -y -t $OMP_NUM_THREADS -x map-ont ../hg38.mmi ../../input/emma/batch_4/$experiment_id"/"$experiment_id"_pass_barcode"$barcode".fastq.gz" > $participant_id"_1_all_reads.sam"

echo \($(date +"%Y-%m-%d %T")\) Sorting and binarising aligned reads
samtools sort -@ $OMP_NUM_THREADS $participant_id"_1_all_reads.sam" > $participant_id"_1_all_reads.bam"

echo \($(date +"%Y-%m-%d %T")\) Indexing sorted reads
samtools index -@ $OMP_NUM_THREADS $participant_id"_1_all_reads.bam"

echo \($(date +"%Y-%m-%d %T")\) Deleting uncompressed reads
rm $participant_id"_1_all_reads.sam"

echo; echo \($(date +"%Y-%m-%d %T")\) done
