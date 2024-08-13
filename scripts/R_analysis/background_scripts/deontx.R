library(Seurat)
library(celda)
library(SingleCellExperiment)
library(ggplot2)

setwd("/data/scratch/projects/punim1441/Project_cortex_organoid_LRsc/analysis/")

directory <- "empty_drops"

################## Define functions ############
#fucntion to run deconx
run_decontX <- function(seurat_obj, background_counts_path) {
  # Extract filtered counts from the Seurat object
  filtered_counts <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))
  
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
  contamination_summary <- summary(sce$decontX_contamination)
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
  ) + labs(title = "old UMAP with Clusters")
  
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
  
  # Display the combined plot
  print(combined_umap)
  
  return(list(seurat_obj = seurat_obj, decontaminated_counts = decontaminated_counts, contamination_summary =contamination_summary))
}

# Function to run decontX over multiple inputs using mapply
run_multiple_decontX <- function(seurat_objs, background_paths, output_names) {
  if (length(seurat_objs) != length(background_paths) || length(seurat_objs) != length(output_names)) {
    stop("All input vectors must have the same length.")
  }
  
  results <- mapply(function(seurat_obj, background_path, output_name) {
    cat("Processing:", output_name, "\n")
    result <- run_decontX(seurat_obj, background_path)
    return(result)
  }, seurat_objs, background_paths, output_names, SIMPLIFY = FALSE)
  
  names(results) <- output_names
  return(results)
}

### fucntion i can use to generate a new dec seuratobject object 
umap_decom_QC <- function(count.matrix, npc = 20, fig_name = '', project = "", MT = 10) {
  
  ###make rst list
  rst_figures <- list()
  rst_table <- data.frame()
  
  #####
  # Initialize Seurat object
  seurat_object <- CreateSeuratObject(counts = count.matrix, project = project)
  
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  # Plot violin plots
  plot(VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")))
  
  vln1 <- VlnPlot(seurat_object, features = c("nFeature_RNA"))
  vln2 <- VlnPlot(seurat_object, features = c("nCount_RNA"))
  vln3 <- VlnPlot(seurat_object, features = c("percent.mt"))
  
  #######
  # Normalize data
  seurat_object <- NormalizeData(seurat_object)
  
  # Identify highly variable features
  seurat_object <- FindVariableFeatures(seurat_object)
  
  # Apply linear transformation
  all_genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all_genes)
  
  # Perform PCA
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  
  # Visualize PCA
  rst_figures <- append(rst_figures, ElbowPlot(seurat_object))
  
  # Cluster cells
  seurat_object <- FindNeighbors(seurat_object, dims = 1:npc)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  
  # Perform UMAP
  seurat_object <- RunUMAP(seurat_object, dims = 1:npc)
  
  # Appending in figure 
  rst_figures <- append(rst_figures, list(DimPlot(seurat_object, reduction = "umap") + labs(color = "Cluster \n(from PCA)", title = '') + theme(text = element_text(size = 10))))
  
  rst_figures <- append(rst_figures, list(
    FeaturePlot(seurat_object, reduction = "umap", features = 'nCount_RNA') + labs(color = "UMI count", title = '') + theme(text = element_text(size = 10)),
    FeaturePlot(seurat_object, reduction = "umap", features = 'nFeature_RNA') + labs(color = str_wrap("Feature count (isoform/gene)", 15), title = '') + theme(text = element_text(size = 10))
  ))
  
  plot_scatter2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering")
  
  plot_pc <- ElbowPlot(seurat_object) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10))
  plot_umap <- grid.arrange(plot_pc,
                            #tableGrob(rst_table),
                            plot_scatter2,
                            rst_figures[[10]], rst_figures[[11]], rst_figures[[12]], vln1, vln2, vln3, ncol = 2, top = textGrob(fig_name))
  
  
  stats_sumary <- rbind("Sample ID" = project,
                        "Cells" = dim(seurat_object)[2],
                        "Median Feature per Cell before filter" = median(seurat_object$nFeature_RNA),
                        "Median Reads per Gene/Isoform" = median(seurat_object$nCount_RNA),
                        "MT Percentage" = MT,
                        "NPCs" = npc,
                        "Median Percent MT before Filter" = median(seurat_object@meta.data[["percent.mt"]]))
  
  tbl_sts2 <- tableGrob(stats_sumary)
  
  write.table(stats_sumary, file = paste0(project, "_stats.csv")) 
  
  grid.newpage()
  grid.draw(tbl_sts2)
  
  list(plot_umap, 
       seurat_object, 
       stats_sumary)
}

# Function to read Seurat objects from a directory
read_seurat_objects <- function(directory) {
  file_paths <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)
  seurat_objs <- lapply(file_paths, readRDS)
  names(seurat_objs) <- gsub("\\.rds$", "", basename(file_paths))
  return(seurat_objs)
}

####### run decontx for all the data ######

# Read Seurat objects from the directory
seurat_objs <- read_seurat_objects(directory)

# Define the paths to the background counts
background_paths <- list.files("background/genes/", 
                               pattern = "geneSymbol",
                               full.names = TRUE)

# Define output names based on Seurat object names
output_names <- names(seurat_objs)

# Run decontX on multiple Seurat objects
results <- run_multiple_decontX(seurat_objs, background_paths, output_names)

saveRDS(results, file = "decontx.rds")




#### so will need see how merging and integration goes here #### 
### perhaps it may be best to just use the docntx file as input into the umap_decom_QC function. 
## will leave this here and reassess ## 

###C1STC
pdf(file = "genes_decon/C1STC_QC.pdf", width = 12, height = 12) 
plotsC1_STC <- plot_umap_decom(decontaminated_counts, npc = 10,fig_name = 'decon_C1_STC (gene counts, Kolf2.1)', project = "C1_STC", MT=10)
dev.off()

decon_C1_STC_umap_object <- plotsC1_STC[[2]]

saveRDS(C1_STC_umap_object, file = "C1_STC_umap_object.rds")

###C4_STC
pdf(file = "genes_decon/C4_STC_QC.pdf", width = 12, height = 12) 
plotsC4_STC <- plot_umap_decom(decontaminated_counts, npc = 10, fig_name = 'decon_C4_STC (gene counts, Kolf2.1)', project = "C4_STC", MT=10)
dev.off()

C4_STC_umap_object <- plotsC4_STC[[2]]
saveRDS(C4_STC_umap_object, file = "C4_STC_umap_object.rds")

###C2Day25
pdf(file = "C2Day25_QC.pdf", width = 12, height = 12) 
plotsC2Day25 <- plot_umap_decom(decontaminated_counts, npc = 10, fig_name = 'C2_Day25 (gene counts, Kolf2.1)', project = "C2_Day25", MT=10)
dev.off()

C2_Day25_umap_object <- plotsC2Day25[[2]]
saveRDS(C2_Day25_umap_object, file = "C2_Day25_umap_object.rds")


plotsC5Day80_decon <- plot_umap_decom(decontaminated_counts, npc = 10,fig_name = 'C5Day80_decon (gene counts, Kolf2.1)', project = "C5_Day80", MT=10)




