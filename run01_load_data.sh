#!/bin/bash

#SBATCH --array=0-61
#SBATCH --job-name=load 
#SBATCH --mail-type=END,FAIL     

#SBATCH --mail-user=yongli.qi@yale.edu
#SBATCH -p day
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH -t 1:00:00

SAMPLE_LIST=($(<SamplesFileNames.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

sample_in=${SAMPLE}

echo "Running on sample  "$SAMPLE >> log

module load R/4.3.0-foss-2022b

Rscript output_seurat_obj.r -i ../00data/$sample_in
