#!/bin/bash -e
#SBATCH -A uoa04084                                               # Project Account
#SBATCH --time=12:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=6                                         # Number of threads (level of paralellism)
#SBATCH --partition=large
#SBATCH --mail-user=ejad042@aucklanduni.ac.nz
#SBATCH --mail-type=END

#SBATCH --output=slurm_output/slurm-%A_%a-IndexingReference.out


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load minimap2/2.27-GCC-12.3.0



cd output

echo \($(date +"%Y-%m-%d %T")\) Indexing reference genome
minimap2 -x map-ont -d hg38.mmi ../input/hg38_test/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz      ## index ref genome (once only)

echo; echo \($(date +"%Y-%m-%d %T")\) done