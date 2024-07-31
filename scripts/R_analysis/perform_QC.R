# load all libraries i need for workflow
library(DoubletFinder)
library(scater)
library(scran)
library(BiocSingular)
library(monocle3)
library(gridExtra)
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(stringr)
library(grid)
library(cowplot)
library(tidyverse)



#setwd("/data/scratch/users/yairp/FLAMES-may1/analysis/seurat_analysis/data/genes/")

###########QC and filtering #######

# get data location for filt CSV objects
org_1A <- read.csv('geneSymbol_org_1A_gene_count.csv', header=T, row.names = 1)
org_1B <- read.csv('geneSymbol_org_1B_gene_count.csv',header=T, row.names = 1)

org_3A <- read.csv('geneSymbol_org_3A_gene_count.csv',header=T, row.names = 1)
org_3B <- read.csv('geneSymbol_org_3B_gene_count.csv',header=T, row.names = 1)
org_3C <- read.csv('geneSymbol_org_3C_gene_count.csv',header=T, row.names = 1)

org_6A <- read.csv('geneSymbol_org_6A_gene_count.csv',header=T, row.names = 1)
org_6B <- read.csv('geneSymbol_org_6B_gene_count.csv',header=T, row.names = 1)
org_6C <- read.csv('geneSymbol_org_6C_gene_count.csv',header=T, row.names = 1)



#setwd("/data/scratch/users/yairp/FLAMES-may1/analysis/seurat_analysis/data/genes_QC/")

### plot UMAPS and sumamry stats from filtered objects 
##UMAP plotting function


#fucntion takes a single cell count matrix -> outputputs sumamry plots and UMAP object and removes doublets.
plot_umap <- function(count.matrix=C1_STC, min.features = 0 , max.features = 100000, max.counts = 10000, min.counts = 1000, npc = 20, cluster_res = 0.7, fig_name = '1', project = "2", MT = 25) {
  
  # Function to calculate min.features and max.features if not provided
  calculate_feature_range <- function(nFeature_RNA) {
    min_feature <- round(mean(nFeature_RNA) - (1.5 * sd(nFeature_RNA)))
    max_feature <- round(mean(nFeature_RNA) + (1.5 * sd(nFeature_RNA)))
    return(list(min_feature = min_feature, max_feature = max_feature))
  }
  
  # If min.features and max.features not provided, calculate them
  if (is.null(min.features) || is.null(max.features)) {
    seurat_obj_org <- CreateSeuratObject(counts = count.matrix, project = project, min.cells = 5, min.features = 1)
    feature_range <- calculate_feature_range(seurat_obj_org$nFeature_RNA)
    min.features <- feature_range$min_feature
    max.features <- feature_range$max_feature
  }
  
  rst_figures <- list()
  rst_table <- data.frame()
  
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
  seurat_object <- CreateSeuratObject(counts = counts_sub, project = project, min.cells = 5, min.features = 1)
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
  #rst_table <- rbind(rst_table, data.frame("Cells" = dim(seurat_object)[2], "Median Feature per Cell" = median(new_seurat_object$nFeature_RNA), "Median Reads per Feature" = median(new_seurat_object$nCount_RNA), row.names = paste0('filter'), check.names = FALSE))
  
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
  
  # Visualize PCA
  rst_figures <- append(rst_figures, ElbowPlot(new_seurat_object))
  
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
  
  # Appending in figure 
  fig_umap <- DimPlot(new_seurat_object, reduction = "umap") + labs(color = "Cluster \n(from PCA)", title = '') + theme(text = element_text(size = 10))
  fig_count <- FeaturePlot(new_seurat_object, reduction = "umap", features = 'nCount_RNA') + labs(color = "UMI count", title = '') + theme(text = element_text(size = 10))
  fig_feature <- FeaturePlot(new_seurat_object, reduction = "umap", features = 'nFeature_RNA') + labs(color = str_wrap("Feature count (isoform/gene)", 15), title = '') + theme(text = element_text(size = 10))
  
  
  plot_scatter2 <- FeatureScatter(new_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering")
  
  
  plot_pc <- ElbowPlot(new_seurat_object) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10))
  plot_umap <- grid.arrange(plot_pc,
                            #tableGrob(rst_table),
                            plot_scatter2, fig_umap, fig_count, fig_feature, vln1, vln2, vln3, ncol = 2, top = textGrob(fig_name))
  
  plot(doublets)
  tbl_sts1 <- tableGrob(statsDoublets)
  grid.newpage()
  grid.draw(tbl_sts1)
  
  stats_sumary <- rbind("Sample ID" = project,
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
  
  write.table(stats_sumary, file = paste0(project, "_stats.csv")) 
  
  grid.newpage()
  grid.draw(tbl_sts2)
  
  list(plot_umap, 
       new_seurat_object, 
       statsDoublets, 
       stats_sumary, 
       new_seurat_object_doublets)
}


#### Generate QC plots and output umap object

###org_1A
pdf(file = "org_1A_QC.pdf", width = 12, height = 12) 
plotsC1_STC <- plot_umap(org_1A, min.features = 2500, max.features = 100000, max.counts =100000, min.counts=2000, npc = 15, cluster_res = 0.5,fig_name = 'org_1A (gene counts, cortex_org)', project = "org_1A", MT=10)
dev.off()

org1A_umap_object <- plotsC1_STC[[2]]
C1_STC_umap_object_doublets <- plotsC1_STC[[5]]

saveRDS(org1A_umap_object, file = "org1A_umap_object.rds")

###org_1B
pdf(file = "org_1B_QC.pdf", width = 12, height = 12) 
plotsorg1_B <- plot_umap(org_1B, min.features = 2500, max.features = 100000, max.counts =100000, min.counts=1000, npc = 15, cluster_res = 0.7,fig_name = 'org_1B (gene counts, cortex_org)', project = "org_1B", MT=10)
dev.off()

plotsorg1_B <- plotsorg1_B[[2]]
saveRDS(plotsorg1_B, file = "org_1B_umap_object.rds")

###org_3A
pdf(file = "org_3A_QC.pdf", width = 12, height = 12) 
plots_org_3A <- plot_umap(org_3A, min.features = 2500, max.features = 10000, max.counts =100000, min.counts=2000,  npc = 15, cluster_res = 0.7,fig_name = 'org_3A (gene counts, cortex_org)', project = "org_3A", MT=10)
dev.off()

plots_org_3A <- plots_org_3A[[2]]
saveRDS(plots_org_3A, file = "org_3A_umap_object.rds")

###org_3B
pdf(file = "org_3B_QC.pdf", width = 12, height = 12) 
plotsorg_3B <- plot_umap(org_3B, min.features = 2500, max.features = 100000, max.counts =100000, min.counts=2000, npc = 15, cluster_res = 0.7,fig_name = 'org_3A (gene counts, cortex_org)', project = "org_3B", MT=10)
dev.off()

plotsorg_3B_umap <- plotsorg_3B[[2]]
saveRDS(plotsorg_3B_umap, file = "org_3B_umap_object.rds")


###org_3C
pdf(file = "org_3C_QC.pdf", width = 12, height = 12) 
plotsC2Day25 <- plot_umap(org_3C, min.features = 2500, max.features = 10000, max.counts =100000, min.counts=2000, npc = 15, cluster_res = 0.7,fig_name = 'org_3C (gene counts, cortex_org)', project = "org_3C", MT=10)
dev.off()

C2_Day25_umap_object <- plotsC2Day25[[2]]
saveRDS(C2_Day25_umap_object, file = "org_3C_umap_object.rds")


###org_6A
pdf(file = "org_6A_QC.pdf", width = 12, height = 12) 
plotsC2Day55 <- plot_umap(org_6A, min.features = 2500, max.features = 10000, max.counts =100000, min.counts=2000, npc = 15, cluster_res = 0.7,fig_name = 'org_6A (gene counts, cortex_org)', project = "org_6A", MT=10)
dev.off()

C2_Day55_umap_object <- plotsC2Day55[[2]]
saveRDS(C2_Day55_umap_object, file = "org_6A_umap_object.rds")

###org_6B
pdf(file = "org_6B_QC.pdf", width = 12, height = 12) 
plotsorg_6B <- plot_umap(org_6B, min.features = 2500, max.features = 10000, max.counts =100000, min.counts=2000, npc = 15, cluster_res = 0.7,fig_name = 'org_6B (gene counts, cortex_org)', project = "org_6B", MT=10)
dev.off()

plotsorg_6B_umap_object <- plotsorg_6B[[2]]
saveRDS(plotsorg_6B_umap_object, file = "org_6B_umap_object.rds")

###org_6C
pdf(file = "org_6C_QC.pdf", width = 12, height = 12) 
plotsorg_6C <- plot_umap(org_6C, min.features = 2500, max.features = 10000, max.counts =100000, min.counts=2000, npc = 15, cluster_res = 0.7,fig_name = 'org_6C (gene counts, cortex_org)', project = "org_6C", MT=10)
dev.off()

plotsorg_6C_umap_object <- plotsorg_6C[[2]]
saveRDS(plotsorg_6C_umap_object, file = "org_6C_umap_object.rds")
