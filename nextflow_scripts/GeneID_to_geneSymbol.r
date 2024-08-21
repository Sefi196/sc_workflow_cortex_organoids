#!/usr/bin/env Rscript

library(dplyr)
library(optparse)

option_list <- list(
    make_option("--data", type="character", help="Path to input data file"),
    make_option("--ref_table", type="character", help="Path to reference table"),
    make_option("--output", type="character", help="Path to output file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Load the data
data_obj <- read.csv(opt$data, header=T)
ref_table <- read.csv(opt$ref_table, header = TRUE)

colnames(data_obj)[1] <- "gene_id"

merged.file = merge(data_obj, ref_table, by.x='gene_id', by.y='gene_id')

duplicates <- duplicated(merged.file$genesymbol)
df_unique <- subset(merged.file, !duplicates)

df_unique <- df_unique[, !colnames(df_unique) %in% "gene_id"]
row.names(df_unique) <- df_unique[, ncol(df_unique)]
df_unique <- df_unique[, -ncol(df_unique)]

# Write output
write.csv(df_unique, opt$output)
