library(FLAMES)
config_file <- "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/resources/config.json"

####Inputs that change depending on sampLe

#####fastq_files 
fastqs = c("/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/fastqs/org_1A.fastq.gz", "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/fastqs/org_3A_fastq.gz") # or fastqs = c("sample1.fq", "sample2.fq")
output = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/outs/batch_1"

####referecne files 
GTF = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/resources/gencode.v41.annotation.gtf"
genome = "/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/resources/true.hg38.analysisSet.fa"
minimap2_dir="/home/yairp/anaconda3/envs/siceLore/bin"


#Run FLAMES multisample pipeline using blaze
sce <- sc_long_multisample_pipeline(fastqs=fastqs, outdir=output, annot=GTF, genome_fa=genome, 
                         #minimap2_dir=minimap2_dir,
                        config_file=config_file, expect_cell_number=c(2000, 2000) )
