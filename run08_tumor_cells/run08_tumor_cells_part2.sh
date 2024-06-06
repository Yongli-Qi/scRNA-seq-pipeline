#!/bin/bash
#SBATCH --job-name=infercnv_HMM3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yongli.qi@yale.edu
#SBATCH -p day
#SBATCH --cpus-per-task=5
#SBATCH --mem=500G
#SBATCH -t 24:00:00

module load R/4.3.0-foss-2022b

export PATH=/home/yq238/miniconda3/bin:$PATH
source activate scRNA

export LD_LIBRARY_PATH=/home/yq238/.conda/envs/scRNA/lib/JAGS/modules-4:$LD_LIBRARY_PATH

cd /home/yq238/project/scRNA_seq/HTAN/08tumor_cells
Rscript infercnv_HMM3.r

