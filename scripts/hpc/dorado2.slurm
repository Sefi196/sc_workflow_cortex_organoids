#!/bin/bash

#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --job-name="org_3B_dorado"
#SBATCH --account="punim1441"
#SBACTH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=yairp@student.unimelb.edu.au
#SBATCH --mail-type=ALL

module load dorado
module load BEDTools/2.30.0
module load pigz

dorado basecaller /data/scratch/projects/punim1441/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 /data/scratch/users/yairp/3B/pod5 --min-qscore 10 --emit-fastq --verbose  > /data/scratch/users/yairp/3B/fastq/org_3B.fastq

#bedtools bamtofastq -i /data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/20240625_Cortex_organoid_scLR_1A_3A/org_1A/20240625_1728_P2S-01556-A_PAW28745_3e141570/fastq/org_1A_calls.bam -fq /data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/20240625_Cortex_organoid_scLR_1A_3A/org_1A/20240625_1728_P2S-01556-A_PAW28745_3e141570/fastq/org_1A.fastq 

pigz /data/scratch/users/yairp/3B/fastq/org_3B.fastq
