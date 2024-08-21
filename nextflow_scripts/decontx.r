suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(celda))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(cowplot))

# Define command-line options
option_list = list(
  make_option(c("--seurat_obj"), type="character", default=NULL, help="Path to Seurat object file (.rds)", metavar="character"),
  make_option(c("--background_counts_path"), type="character", default=NULL, help="Path to background counts file (.csv)", metavar="character"),
  make_option(c("--sample_id"), type="character", default=NULL, help="Output file name for the results (.rds)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$seurat_obj) || is.null(opt$background_counts_path) || is.null(opt$sample_id)) {
  stop("Please provide all arguments: --seurat_obj, --background_counts_path, and --sample_id")
}

# Function to run decontX on a single Seurat object
run_decontX <- function(seurat_obj_path, background_counts_path, sample_id) {
  # Load the Seurat object
  seurat_obj <- readRDS(seurat_obj_path)
  
  filtered_counts <- as.matrix(GetAssayData(seurat_obj, layer = "counts"))
  
  # Read background counts
  raw_counts <- as.matrix(read.csv(background_counts_path, row.names = 1))
  
  # Get cluster info from Seurat object
  cluster_info <- setNames(seurat_obj$seurat_clusters, colnames(seurat_obj))
  
  # Find common genes
  common_genes <- intersect(rownames(filtered_counts), rownames(raw_counts))
  raw_counts <- raw_counts[common_genes, ]
  filtered_counts <- filtered_counts[common_genes, ]
  
  # Create SingleCellExperiment objects
  sce_raw <- SingleCellExperiment(list(counts = raw_counts))
  sce_object <- SingleCellExperiment(list(counts = filtered_counts))
  
  # Run decontX with background
  sce <- decontX(sce_object, z = cluster_info, background = sce_raw)
  
  # Summarize contamination levels
  contamination_summary <- as.array(summary(sce$decontX_contamination))
  print(contamination_summary)
  
  # Add contamination levels to Seurat object metadata
  contamination <- colData(sce)$decontX_contamination
  seurat_obj <- AddMetaData(seurat_obj, metadata = contamination, col.name = "decontX_contamination")
  
  # Extract decontaminated counts from SCE object
  decontaminated_counts <- assay(sce, "decontXcounts")
  decontaminated_counts <- as.matrix(decontaminated_counts)
  
  # Create a new assay with decontaminated counts and add it to Seurat object
  new_assay <- CreateAssayObject(counts = decontaminated_counts)
  seurat_obj[["decontaminated"]] <- new_assay
  
  clusters_umap_orig <- DimPlot(
    object = seurat_obj,
    group.by = "seurat_clusters",
    reduction = "umap",
    label = TRUE,
    pt.size = 0.5
  ) + labs(title = "Old UMAP with Clusters")
  
  # Plot UMAP with contamination levels
  contamination_umap <- FeaturePlot(
    object = seurat_obj, 
    features = "decontX_contamination", 
    reduction = "umap"
  ) + labs(title = "UMAP Colored by decontX_contamination")
  
  DefaultAssay(seurat_obj) <- "decontaminated"
  
  # Normalization, variable feature selection, and scaling
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  # PCA and clustering
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.7)
  
  # UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  
  # Plot UMAP with updated clusters
  clusters_umap <- DimPlot(
    object = seurat_obj,
    group.by = "seurat_clusters",
    reduction = "umap",
    label = TRUE,
    pt.size = 0.5
  ) + labs(title = "UMAP with Corrected Clusters")
  
  # Combine plots
  combined_umap <- cowplot::plot_grid(clusters_umap_orig, contamination_umap, clusters_umap, ncol = 3)
  
  cat("Making plots\n")
  # Save the combined plot as a PDF
  pdf(file = paste0(sample_id, "_decontx_plots.pdf"), width = 18, height = 6)
  print(combined_umap)
  dev.off()

  cat("Saving seurat obj\n")
  # Save the Seurat object
  saveRDS(seurat_obj, file = paste0(sample_id, "_decontx_seurat_obj.rds"))

  # Save decontaminated counts and contamination summary
  cat("Saving decontx counts\n")
  write.csv(decontaminated_counts, paste0(sample_id, "_decontx_counts.csv"))
  
  # Print a message indicating that the contamination summary is being saved
  cat("Saving contamination summary\n")
  
  # Ensure contamination_summary is a data frame
  contamination_summary_df <- as.data.frame(contamination_summary)
  write.table(contamination_summary_df, file = paste0(sample_id, "_contamination_summary.txt"))

  # Optionally return the results
  #return(list(seurat_obj = seurat_obj, decontaminated_counts = decontaminated_counts, contamination_summary = contamination_summary))
}

# Run decontX on the input Seurat object and background counts file
result <- run_decontX(opt$seurat_obj, opt$background_counts_path, opt$sample_id)

cat("Results saved to:", paste0(opt$sample_id, "_decontx_seurat_obj.rds"), "\n")
