#### isopod testing on cortex diffs ####


library(isopod)
library(data.table)
library(plyr)
library("dplyr")
library(UpSetR)
library(tidyr)
library(tidyverse)

count_matrix <- as.data.frame(GetAssayData(obj, assay = "iso", layer = "counts"))
# Convert matrix to data frame and add row names as a new column
count_df <- cbind(transcript_id = rownames(count_matrix), count_matrix)
# Reset row names
rownames(count_df) <- NULL


### need to add back in the gene_id col. 
ENST_to_ENSG <- read.csv("gene_transcript_isd.csv")

count_table <- merge(count_df, ENST_to_ENSG, by= "transcript_id")

counts_table <-read.csv("/Users/yairp/Library/CloudStorage/Dropbox/My Mac (5160L-152961-M)/Documents/Projects/Single_Cell_10x/BLAZE/Spartan_working_dir/script/analysis/data/trans_count_B+F_matt.csv")
cell_groups <- obj$



permutation_results <- run_everything(counts_table, 
                                      cell_groups, 
                                      transcript_id_colname = 'transcript_id', 
                                      gene_id_colname = 'gene_id',
                                      cell_labels_colname = 'groups', 
                                      cell_group_to_analyse = 'A549', 
                                      output_folder = 'permutation_results')



