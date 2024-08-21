suppressPackageStartupMessages({
library(patchwork)
library(Matrix)
library(ggplot2)
library(stringr)
library(grid)
library(cowplot)
library(tidyverse)
library(optparse)
library(Seurat)
library(DoubletFinder)
library(gridExtra)
})

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input CSV file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, help="Sample ID", metavar="character"),
  make_option(c("-p", "--npc"), type="integer", default=20, help="Number of principal components", metavar="integer"),
  make_option(c("-r", "--cluster_res"), type="numeric", default=0.7, help="Clustering resolution", metavar="numeric"),
  make_option(c("-m", "--mt"), type="numeric", default=10, help="Mitochondrial percentage threshold", metavar="numeric"),
  make_option(c("--min_counts"), type="integer", default=1000, help="Minimum counts", metavar="integer"),
  make_option(c("--max_counts"), type="integer", default=100000, help="Maximum counts", metavar="integer"),
  make_option(c("--min_features"), type="integer", default=NULL, help="Minimum features", metavar="integer"),
  make_option(c("--max_features"), type="integer", default=NULL, help="Maximum features", metavar="integer")
)

opt <- parse_args(OptionParser(option_list = option_list))


###########QC and filtering #######
### plot UMAPS and sumamry stats from filtered objects 
##UMAP plotting fucntion
#fucntion takes a single cell count matrix -> outputputs sumamry plots and UMAP object and removes doublets.
plot_umap <- function(count.matrix=C1_STC, min.features = NULL, max.features = NULL, max.counts = 10000, min.counts = 10000, npc = 20, cluster_res = 0.7, sample = '1', MT = 10) {
  
  # Function to calculate min.features and max.features if not provided
  calculate_feature_range <- function(nFeature_RNA) {
    min_feature <- round(mean(nFeature_RNA) - (1.5 * sd(nFeature_RNA)))
    max_feature <- round(mean(nFeature_RNA) + (1.5 * sd(nFeature_RNA)))
    return(list(min_feature = min_feature, max_feature = max_feature))
  }
  
  # If min.features and max.features not provided, calculate them
  if (is.null(min.features) || is.null(max.features)) {
    seurat_obj_org <- CreateSeuratObject(counts = count.matrix, project = sample, min.cells = 5, min.features = 1)
    feature_range <- calculate_feature_range(seurat_obj_org$nFeature_RNA)
    min.features <- feature_range$min_feature
    max.features <- feature_range$max_feature
  }
  
  pdf(file = paste0(sample, "_plots.pdf"), width = 12, height = 12) 

 
  # Calculate the percentage of cells expressing each gene
  gene_percent_expression <- rowMeans(count.matrix > 0) * 100
  
  # Select genes expressed in at least 1% of cells
  genes_filter <- names(gene_percent_expression[gene_percent_expression > 1])
  
  # Filter counts
  counts_sub <- count.matrix[genes_filter, ]
  
  # Record the number of features removed
  removed_features <- dim(count.matrix)[1] - length(genes_filter)
  
  #####
  # Initialize Seurat object
  seurat_object <- CreateSeuratObject(counts = counts_sub, project = sample, min.cells = 5, min.features = 1)
  #rst_table <- rbind(rst_table, data.frame("Cells" = dim(seurat_object)[2], "Median Feature per Cell" = median(seurat_object$nFeature_RNA), "Median Reads per Feature" = median(seurat_object$nCount_RNA), row.names = paste0('No filter'), check.names = FALSE))
  plot_scatter1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell BEFORE filtering")
  
  plot(plot_scatter1)
  
  seurat_object[["joined"]] <- JoinLayers(seurat_object[["RNA"]])
  
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  # Plot violin plots
  plot(VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")))
  
  # Remove unwanted cells
  new_seurat_object <- subset(seurat_object, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < MT & nCount_RNA < max.counts & nCount_RNA > min.counts) 
  
  # Plot violin plots again after filtering
  plot(VlnPlot(new_seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")))
  
  vln1 <- VlnPlot(new_seurat_object, features = c("nFeature_RNA"))
  vln2 <- VlnPlot(new_seurat_object, features = c("nCount_RNA"))
  vln3 <- VlnPlot(new_seurat_object, features = c("percent.mt"))
  
  #######
  # Normalize data
  new_seurat_object <- NormalizeData(new_seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identify highly variable features
  new_seurat_object <- FindVariableFeatures(new_seurat_object, selection.method = "vst", nfeatures = 2000)
  
  # Apply linear transformation
  all_genes <- rownames(new_seurat_object)
  new_seurat_object <- ScaleData(new_seurat_object, features = all_genes)
  
  # Perform PCA
  new_seurat_object <- RunPCA(new_seurat_object, features = VariableFeatures(object = new_seurat_object))
  
  # Cluster cells
  new_seurat_object <- FindNeighbors(new_seurat_object, dims = 1:npc)
  new_seurat_object <- FindClusters(new_seurat_object, resolution = cluster_res)
  
  # Perform UMAP
  new_seurat_object <- RunUMAP(new_seurat_object, dims = 1:npc)
  
  ### Filter out doublets (remember to modify doublet rate if samples have variable target cells)
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_pbmc <- paramSweep(new_seurat_object, PCs = 1:20, sct = FALSE)
  sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
  bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
  
  pK <- bcmvn_pbmc %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- new_seurat_object@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.039 * nrow(new_seurat_object@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run doubletFinder 
  new_seurat_object <- doubletFinder(new_seurat_object, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  colnames(new_seurat_object@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(new_seurat_object@meta.data))
  
  # Summary doublets
  statsDoublets <- new_seurat_object@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(Median_nCount_RNA = median(nCount_RNA), Median_nFeature_RNA = median(nFeature_RNA), Count = n())
  
  # Visualize doublets
  doublets <- DimPlot(new_seurat_object, reduction = 'umap', group.by = "DF.classifications")
  
  ### i want to save the seurat object with doublets listed 
  new_seurat_object_doublets <- new_seurat_object
  
  new_seurat_object <- subset(new_seurat_object, subset = DF.classifications == 'Singlet')
  
  # figures
  umap_plot <- DimPlot(new_seurat_object, reduction = "umap") + labs(color = "Cluster \n(from PCA)", title = '') + theme(text = element_text(size = 10))
  
  Fplot1 <- FeaturePlot(new_seurat_object, reduction = "umap", features = 'nCount_RNA') + labs(color = "UMI count", title = '') + theme(text = element_text(size = 10))
  Fplot2 <- FeaturePlot(new_seurat_object, reduction = "umap", features = 'nFeature_RNA') + labs(color = str_wrap("Feature count (gene)", 15), title = '') + theme(text = element_text(size = 10))
  
  plot_scatter2 <- FeatureScatter(new_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering")
  
  plot_pc <- ElbowPlot(new_seurat_object) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10))
  
  ggplot_list <- list(plot_pc, plot_scatter2, umap_plot, Fplot1, Fplot2, vln1, vln2, vln3)
  combined_plots <- plot_grid(plotlist = ggplot_list, ncol = 2)

  plot(combined_plots)
  plot(doublets)
  
  tbl_sts1 <- tableGrob(statsDoublets)
  grid.newpage()
  grid.draw(tbl_sts1)
  
  stats_sumary <- rbind("Sample ID" = sample,
                        "Cells_before_filter" = dim(seurat_object)[2],
                        "Cells_after_filter" = dim(new_seurat_object)[2],
                        "Median Feature per Cell before filter" = median(seurat_object$nFeature_RNA),
                        "Median Reads per Gene/Isoform before filter" = median(seurat_object$nCount_RNA),
                        "Median Feature per Cell" = median(new_seurat_object$nFeature_RNA),
                        "Median Reads per Gene/Isoform" = median(new_seurat_object$nCount_RNA),
                        "Max Features" = max.features,
                        "Min Features" = min.features,
                        "Min Counts" = min.counts,
                        "Max Counts" = max.counts,
                        "MT Percentage" = MT,
                        "NPCs" = npc,
                        "Median Percent MT before Filter" = median(seurat_object@meta.data[["percent.mt"]]),
                        "Median Percent MT after Filter" = median(new_seurat_object@meta.data[["percent.mt"]]),
                        "Removed Features" = removed_features)
  
  tbl_sts2 <- tableGrob(stats_sumary)
  
  grid.newpage()
  grid.draw(tbl_sts2)
  
  dev.off()
  
  
  cat("saving seurat objects\n")

  saveRDS(new_seurat_object, file = paste0(sample, "_umap_object.rds"))
  saveRDS(new_seurat_object_doublets, file = paste0(sample, "_with_doublets_umap_object.rds"))
  write.table(stats_sumary, file = paste0(sample, "_stats.csv")) 
}


# Main function
main <- function() {
  if (is.null(opt$input)) {
    stop("Input file is required")
  }
  
  # Load data
  count_matrix <- read.csv(opt$input, header=T, row.names = 1)
  
  # Run UMAP plot function
  result <- plot_umap(count.matrix = count_matrix,
                      min.features = opt$min_features,
                      max.features = opt$max_features,
                      max.counts = opt$max_counts,
                      min.counts = opt$min_counts,
                      npc = opt$npc,
                      cluster_res = opt$cluster_res,
                      sample = opt$sample,
                      MT = opt$mt)
}

# Run the main function
main()