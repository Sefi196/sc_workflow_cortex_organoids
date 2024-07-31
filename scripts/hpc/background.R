library(FLAMES)
config_file <- "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/resources/background.json"

####Inputs that change depending on sampLe

#####fastq_files 
fastqs = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/fastqs/"
output = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/outs/background_batch1_2"
barcodes_file = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/outs/background_batch1_2/backlist"


GTF = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/resources/gencode.v41.annotation.gtf"
genome = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/resources/true.hg38.analysisSet.fa"
minimap2_dir="/home/yairp/anaconda3/envs/siceLore/bin"


#Run FLAMES multisample pipeline using blaze
sce <- sc_long_multisample_pipeline(fastqs=fastqs, outdir=output, annot=GTF, genome_fa=genome, barcodes_file = barcodes_file, config_file=config_file, expect_cell_number=c(2000, 2000) )
