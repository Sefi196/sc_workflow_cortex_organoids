#####classifaction
# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

####
# define functions
perform_sctype_analysis <- function(seurat_obj, db_, tissue, gs_removal_list = c(), 
                                    metadat_col_prefix = "db_prefix", figure_prefix = "fig_name",
                                    cluster_res = "RNA_snn_res.0.5", output_file = "") {
  # Prepare gene sets
  gs_list <- gene_sets_prepare(db_, tissue)
  
  # Remove specified gene sets
  for (gs in gs_removal_list) {
    gs_list[["gs_positive"]][[gs]] <- NULL
  }
  
  # Calculate sctype scores
  es.max <- sctype_score(scRNAseqData = seurat_obj@assays$RNA$scale.data, scaled = TRUE, 
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # Set identities in Seurat object
  Idents(seurat_obj) <- cluster_res
  
  # Merge by cluster
  cL_results <- do.call("rbind", lapply(unique(seurat_obj@meta.data[[cluster_res]]), function(cl) {
    es.max.cl <- sort(rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_res]] == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data[[cluster_res]] == cl)), 10)
  }))
  
  sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  # Set low-confident clusters to "Unknown"
  sctype_scores$scores <- as.numeric(sctype_scores$scores)
  sctype_scores$type[sctype_scores$scores < sctype_scores$ncells / 4] <- "Unknown"
  print(sctype_scores[, 1:3])
  
  # Overlay the labels
  seurat_obj@meta.data[[metadat_col_prefix]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j,]
    seurat_obj@meta.data[[metadat_col_prefix]][seurat_obj@meta.data[[cluster_res]] == j] <- as.character(cl_type$type[1])
  }
  
  # Plotting
  pclass <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = metadat_col_prefix)
  print(pclass)
  
  # Save the plot to a PDF
  pdf(file = paste0(figure_prefix, "_", metadat_col_prefix, "_sctype_genes.pdf"), width = 8, height = 8)
  print(pclass)
  dev.off()
  
  # Save the updated Seurat object to an RDS file
  if (output_file != "") {
    saveRDS(seurat_obj, file = paste0(output_file, ".rds"))
  }
  
  # Return the updated Seurat object
  return(seurat_obj)
}
####
# Define variable

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
#db_ <- "/data/gpfs/projects/punim1441/FLAMES_202311/resources/6mOrgMarkers.xlsx"
tissue <- "Brain"
gs_removal_list <- c("Tanycytes")

org1A_umap_object <- perform_sctype_analysis(org1A_umap_object, db_, tissue, gs_removal_list, 
                        metadat_col_prefix ="manveer", figure_prefix = "org_1A",
                        output_file = "org1A_umap_object", cluster_res = "RNA_snn_res.0.5")


org_1B_umap_object <- perform_sctype_analysis(org_1B_umap_object, db_, tissue, gs_removal_list, 
                                              metadat_col_prefix ="scytpe_db", figure_prefix = "org_1B",
                                              output_file = "org1B_umap_object", cluster_res = "RNA_snn_res.0.7")

org_3A_umap_object <- perform_sctype_analysis(org_3A_umap_object, db_, tissue, gs_removal_list, 
                                             metadat_col_prefix ="scytpe_db", figure_prefix = "org_3A",
                                             output_file = "org3A_umap_object", cluster_res = "RNA_snn_res.0.7")

org_3B_umap_object <- perform_sctype_analysis(org_3B_umap_object, db_, tissue, gs_removal_list, 
                                              metadat_col_prefix ="scytpe_db", figure_prefix = "org_3B",
                                              output_file = "org3B_umap_object", cluster_res = "RNA_snn_res.0.7")

org_3C_umap_object <- perform_sctype_analysis(org_3C_umap_object, db_, tissue, gs_removal_list, 
                                              metadat_col_prefix ="scytpe_db", figure_prefix = "org_3C",
                                              output_file = "org3C_umap_object", cluster_res = "RNA_snn_res.0.7")

org_6A_umap_object <- perform_sctype_analysis(org6A_umap_object, db_, tissue, gs_removal_list, 
                                              metadat_col_prefix ="scytpe_db", figure_prefix = "org_6A",
                                              output_file = "org6A_umap_object", cluster_res = "RNA_snn_res.0.7")

org_6B_umap_object <- perform_sctype_analysis(org_6B_umap_object, db_, tissue, gs_removal_list, 
                                              metadat_col_prefix ="scytpe_db", figure_prefix = "org_6B",
                                              output_file = "org6B_umap_object", cluster_res = "RNA_snn_res.0.7")

org_6C_umap_object <- perform_sctype_analysis(org_6C_umap_object, db_, tissue, gs_removal_list, 
                                              metadat_col_prefix ="scytpe_db", figure_prefix = "org_6C",
                                              output_file = "org6C_umap_object", cluster_res = "RNA_snn_res.0.7")





