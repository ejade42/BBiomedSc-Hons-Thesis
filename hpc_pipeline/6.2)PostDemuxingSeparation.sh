#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=12:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END
#SBATCH --array=0-8
#SBATCH --output=slurm_output/slurm-%A_%a-PostDemuxingSeparation.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load SAMtools/1.19-GCC-12.3.0
module load minimap2/2.27-GCC-12.3.0


i=$SLURM_ARRAY_TASK_ID

cd dorado


if [ $i -eq "0" ]; then
    ## Deal with 20JK4316D (split across PGAXHC230412 (some) and PGAXHX240013 (all))
    participant_id="20JK4316D"
    mkdir -p $participant_id

    experiment_id="PGAXHX240013"
    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_1.bam" $experiment_id"/"$experiment_id"_2_calls.bam"

    experiment_id="PGAXHC230412"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_2.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode17.bam"

    samtools merge -f -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $participant_id"/"$participant_id"_1_calls_exp_1.bam" $participant_id"/"$participant_id"_1_calls_exp_2.bam"
fi



if [ $i -eq "1" ]; then
    ## Deal with 21LE9817N (PGAXHC240012)
    participant_id="21LE9817N"
    experiment_id="PGAXHC240012"
    mkdir -p $participant_id

    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $experiment_id"/"$experiment_id"_2_calls.bam"
fi



if [ $i -eq "2" ]; then
    ## Deal with 22NZ3355R (PGAXHC230412)
    participant_id="22NZ3355R"
    experiment_id="PGAXHC230412"
    mkdir -p $participant_id

    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode19.bam"
fi



if [ $i -eq "3" ]; then
    ## Deal with 311 (PGAXOW240379/PGAXOW240412 barcode NB14)
    participant_id="311"
    mkdir -p $participant_id
    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id
    
    experiment_id="PGAXOW240379"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_1.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode14.bam"
    
    experiment_id="PGAXOW240412"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_2.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode14.bam"
    
    samtools merge -f -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $participant_id"/"$participant_id"_1_calls_exp_1.bam" $participant_id"/"$participant_id"_1_calls_exp_2.bam"
fi



if [ $i -eq "4" ]; then
    ## Deal with 313 (PGAXOW240380/PGAXOW240413 barcode NB11)
    participant_id="313"
    mkdir -p $participant_id
    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id
    
    experiment_id="PGAXOW240380"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_1.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode11.bam"
    
    experiment_id="PGAXOW240413"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_2.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode11.bam"
    
    samtools merge -f -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $participant_id"/"$participant_id"_1_calls_exp_1.bam" $participant_id"/"$participant_id"_1_calls_exp_2.bam"
fi



if [ $i -eq "5" ]; then
    ## Deal with 317 (PGAXOW240380/PGAXOW240413 barcode NB12)
    participant_id="317"
    mkdir -p $participant_id
    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id

    experiment_id="PGAXOW240380"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_1.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode12.bam"
    
    experiment_id="PGAXOW240413"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_2.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode12.bam"
    
    samtools merge -f -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $participant_id"/"$participant_id"_1_calls_exp_1.bam" $participant_id"/"$participant_id"_1_calls_exp_2.bam"
fi



if [ $i -eq "6" ]; then
    ## Deal with 318 (PGAXOW240379/PGAXOW240412 barcode NB10)
    participant_id="318"
    mkdir -p $participant_id
    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id

    experiment_id="PGAXOW240379"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_1.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode10.bam"
    
    experiment_id="PGAXOW240412"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_2.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode10.bam"
    
    samtools merge -f -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $participant_id"/"$participant_id"_1_calls_exp_1.bam" $participant_id"/"$participant_id"_1_calls_exp_2.bam"
fi



if [ $i -eq "7" ]; then
    ## Deal with 319 (PGAXOW240380/PGAXOW240413 barcode NB13)
    participant_id="319"
    mkdir -p $participant_id
    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id
    
    experiment_id="PGAXOW240380"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_1.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode13.bam"
    
    experiment_id="PGAXOW240413"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_2.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode13.bam"
    
    samtools merge -f -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $participant_id"/"$participant_id"_1_calls_exp_1.bam" $participant_id"/"$participant_id"_1_calls_exp_2.bam"
fi



if [ $i -eq "8" ]; then
    ## Deal with 327 (PGAXOW240379/PGAXOW240412 barcode NB09)
    participant_id="327"
    mkdir -p $participant_id
    echo \($(date +"%Y-%m-%d %T")\) extracting reads for participant $participant_id
    
    experiment_id="PGAXOW240379"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_1.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode09.bam"
    
    experiment_id="PGAXOW240412"
    samtools view -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_1_calls_exp_2.bam" $experiment_id"/"$experiment_id"_demuxed/SQK-NBD114-24_barcode09.bam"
    
    samtools merge -f -@ $OMP_NUM_THREADS -o $participant_id"/"$participant_id"_2_calls.bam" $participant_id"/"$participant_id"_1_calls_exp_1.bam" $participant_id"/"$participant_id"_1_calls_exp_2.bam"
fi



echo \($(date +"%Y-%m-%d %T")\) done
