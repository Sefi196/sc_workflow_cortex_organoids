#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --job-name="all_togehter"
#SBATCH --account="punim1441"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=800GB
#SBATCH --mail-user=yairp@student.unimelb.edu.au
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --time=7-0:0:00

module load R/4.4.0
module load minimap2/2.24
module load SAMtools/1.13

Rscript /data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/outs/all_together/FLAMES_all_together.R
