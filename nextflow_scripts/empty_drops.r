# Load required libraries
suppressWarnings({
    library(Seurat)
    library(DropletUtils)
    library(dplyr)
    library(ggplot2)
    library(gridExtra)
    library(data.table)
    library(tibble)
    library(BiocParallel)
    library(grid) 
    })

# Main script to parse command-line arguments and run the function
args <- commandArgs(trailingOnly = TRUE)
output_path <- args[1]
gene_count_file <- args[2]
empty_drops_file <- args[3]
output_seurat_file <- args[4]
fdr_threshold <- as.numeric(args[5])
lower <- as.numeric(args[6])


# Define the function to perform empty drops analysis
perform_empty_drops_analysis <- function(output_path, gene_count_file, empty_drops_file, output_seurat_file, fdr_threshold = 0.001, lower = 100) {
  
  # Create the output directory if it does not exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # Read in data
  print("reading in data")
  df <- read.csv(gene_count_file, row.names = 1)
  df_emptydrops <- read.csv(empty_drops_file, row.names = 1)
  
  # Combine the dataframes by row names
  combined_df <- merge(df, df_emptydrops, by = "row.names", all = TRUE)
  rownames(combined_df) <- combined_df[, 1]
  combined_df[, 1] <- NULL
  combined_df[is.na(combined_df)] <- 0
  
  # Perform standard pre-processing before empty drops analysis
  seurat_obj <- CreateSeuratObject(counts = df, project = "pbmc3k", min.cells = 3, min.features = 20)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 20 & nFeature_RNA < 100000 & percent.mt < 100)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  DimPlot(seurat_obj, reduction = "umap")
  
  # Define function to make dgCMatrix from combined counts
  makedgcmatrix <- function(count.matrix) {
    seurat_object <- CreateSeuratObject(counts = count.matrix, project = "singlecell", min.cells = 3, min.features = 0.5)
    list(seurat_object[["RNA"]]$counts)
  }
  
  # Empty Drops Analysis
  print("Make dgCMatrix")
  outs.ddcmatrix <- makedgcmatrix(combined_df)[[1]]
  br.out <- barcodeRanks(outs.ddcmatrix)
  
  e.out <- emptyDrops(outs.ddcmatrix, lower = lower, niters = 10000, test.ambient = TRUE, BPPARAM = SerialParam())
  is.cell <- e.out$FDR < fdr_threshold
  
  # Create a dataframe with FDR of TRUE cells
  is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
  is.true.cell_CR <- is.true.cell_CR %>% filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
  is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")
  
  # Function for retrieving the Seurat cells and cluster in dataframe
  overlap_true_cell <- function(seurat_object) {
    seurat_cluster.df <- as.data.frame(seurat_object$seurat_clusters)
    seurat_cluster.df <- tibble::rownames_to_column(seurat_cluster.df, "cell_id")
    seurat_cluster.df
  }
  
  # Obtain cluster dataframe from Seurat object
  overlap_CR <- overlap_true_cell(seurat_obj)
  
  # Check overlaps between Seurat object and true cells
  summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)
  
  # Function to add metadata to Seurat object
  print("Finding empty drops and adding to metadata")
  True.cells <- function(e.out) {
    cells <- as.data.frame(e.out@rownames)
    fdr <- as.data.frame(e.out$FDR)
    T.F.cells <- cbind(cells, fdr)
    T.F.cells <- data.frame(T.F.cells[,-1], row.names = T.F.cells[,1])
    setnames(T.F.cells, c('FDR'))
    T.F.cells %>%
      mutate(FDR = case_when(FDR < fdr_threshold ~ "Cells", FDR > fdr_threshold ~ "Empty_drops"))
  }
  
  cells_CR <- True.cells(e.out)
  seurat_obj <- AddMetaData(seurat_obj, metadata = cells_CR, col.name = 'is.cell')
  
  # Create a ggplot object
  rankplot <- ggplot(br.out, aes(x = rank, y = total)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total") +
    geom_line(aes(y = fitted), color = "red", linetype = "solid") +
    geom_hline(yintercept = metadata(br.out)$knee, color = "dodgerblue", linetype = "dashed") +
    geom_hline(yintercept = metadata(br.out)$inflection, color = "forestgreen", linetype = "dashed") +
    theme(
      legend.position = "bottomleft"
    ) +
    guides(colour = guide_legend(override.aes = list(linetype = c("dashed", "dashed")))) +
    annotate("text", x = Inf, y = metadata(br.out)$knee, label = "knee", color = "dodgerblue", vjust = -1, hjust = 1) +
    annotate("text", x = Inf, y = metadata(br.out)$inflection, label = "inflection", color = "forestgreen", vjust = -1, hjust = 1)
  
  # Summary table
  cell_counts <- as.data.frame(table(seurat_obj@meta.data$is.cell))
  count_true_cells <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0, cell_counts$Freq[cell_counts$Var1 == "Cells"], 0)
  count_empty_drops <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0, cell_counts$Freq[cell_counts$Var1 == "Empty_drops"], 0)
  summary_table <- data.frame(
    Description = c('fdr', 'lower Counts', 'number of true cells', 'number of empty drops'),
    Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops)
  )
  summary_grob <- tableGrob(summary_table, rows = NULL, cols = NULL)
  
  # Create the combined plot
  plot1 <- grid.arrange(
    rankplot,
    DimPlot(seurat_obj, reduction = "umap", group.by = 'is.cell') + 
      labs(color = "is.cell", title = 'Seurat Object') + 
      theme(text = element_text(size = 10), plot.background = element_rect(fill = "white")),
    FeaturePlot(seurat_obj, features = "nCount_RNA") + 
      theme(plot.background = element_rect(fill = "white")),
    FeaturePlot(seurat_obj, features = "nFeature_RNA") + 
      theme(plot.background = element_rect(fill = "white")),
    summary_grob,
    ncol = 2,
    top = textGrob('Empty drops vs real cells')
  )
  
  print("Saving data")
  # Output the plot and summary stats
  pdf(file = file.path(output_path, paste0(output_seurat_file, "_plots.pdf")), width = 6, height = 6, bg = "white")
  plot(plot1)
  dev.off()
  
  # Subset the Seurat object to remove cells marked as empty drops
  seurat_obj_rm_empty <- subset(seurat_obj, subset = is.cell == 'Cells')
  
  # Save the Seurat objects
  saveRDS(seurat_obj, file = file.path(output_path, paste0("with_empty_", output_seurat_file, ".rds")))
  saveRDS(seurat_obj_rm_empty, file = file.path(output_path, paste0("removed_empty_", output_seurat_file, ".rds")))
}


#run fucntion
perform_empty_drops_analysis(output_path, gene_count_file, empty_drops_file, output_seurat_file, fdr_threshold, lower)