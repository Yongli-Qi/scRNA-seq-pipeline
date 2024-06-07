###batch script to submit the job###

#!/bin/bash

#SBATCH --array=0-51
#SBATCH --job-name=doubletdetection
#SBATCH --mail-type=END,FAIL     

#SBATCH --mail-user=yongli.qi@yale.edu
#SBATCH -p day
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH -t 1:00:00

SAMPLE_LIST=($(<SamplesFileNames.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

sample_in=/home/yq238/project/scRNA_seq/HTAN/00data/${SAMPLE}
sample_in2=/home/yq238/project/scRNA_seq/HTAN/02seurat_objs_filtered/${SAMPLE}_cells_after_filtering.txt
sample_out=${SAMPLE}_non_doublets_cells.txt
sample_out2=${SAMPLE}_doublets_scores.txt

export PATH=/home/yq238/miniconda3/bin:$PATH
echo "Running on sample  "$SAMPLE >> log
python doubletdetection_run.py $sample_in $sample_in2 $sample_out $sample_out2
