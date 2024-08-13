library(dplyr)
library(Seurat)
library(tidyr)
library(gridExtra)

#read in mapping file 
mapping_file <- read.csv("/data/projects/punim1441/FLAMES_202311/resources/mapping_file_find_gene_transcript.csv")
### define fucntion for modifying coutnmmatrix
## this fucntion cats together gene_id and gene symbol so each row of the count matrix is names ESTID_genesymbol 
modify_row_names <- function(count_matrix, mapping_df) {
  # Ensure the input is a data frame
  count_matrix <- read.csv(count_matrix, header = TRUE, row.names = 1)
  
  # Modify the row names by concatenating transcript_id and gene_id
  count_matrix <- cbind(transcript_id = rownames(count_matrix), count_matrix)
  rownames(count_matrix) <- NULL
  
  # Replace ENSG_IDs with gene symbols
  count_matrix <- merge(count_matrix, mapping_df, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
  
  # Update row names to include gene symbol instead of ENSG_ID
  row.names(count_matrix) <- paste0(count_matrix$transcript_id, "_", count_matrix$GeneSymbol)
  
  # Remove original transcript_id and gene_id columns
  count_matrix$transcript_id <- NULL
  
  # Filter out rows where gene_symbol contains "BambuGene"
  count_matrix <- count_matrix[!grepl("BambuGene", count_matrix$gene_id), ]
  count_matrix$gene_id <- NULL
  
  
  # Return the modified data frame
  return(count_matrix)
}




A <- modify_row_names("org_1A_matched_reads_dedup_transcript_count.csv", mapping_file)
B <- modify_row_names('org_1B_matched_reads_dedup_transcript_count.csv', mapping_file)
C <- modify_row_names('org_3A_matched_reads_dedup_transcript_count.csv',mapping_file)
D <- modify_row_names('org_3B_matched_reads_dedup_transcript_count.csv', mapping_file)
E <- modify_row_names('org_3C_matched_reads_dedup_transcript_count.csv', mapping_file)
FF <- modify_row_names('org_6A_matched_reads_dedup_transcript_count.csv', mapping_file)
G <- modify_row_names('org_6B_matched_reads_dedup_transcript_count.csv', mapping_file)
H <- modify_row_names('org_6C_matched_reads_dedup_transcript_count.csv', mapping_file)

#Make serat objects not filterd
org_1A_umap_object_iso <- CreateSeuratObject(counts = A, project = "org_1A")
org_1B_umap_object_iso <- CreateSeuratObject(counts = B, project = "org_1B")
org_3A_umap_object_iso <- CreateSeuratObject(counts = C , project = "org_3A")
org_3B_umap_object_iso <- CreateSeuratObject(counts = D, project = "org_3B")
org_3C_umap_object_iso <- CreateSeuratObject(counts = E, project = "org_3C")
org_6A_umap_object_iso <- CreateSeuratObject(counts = FF , project = "org_6A")
org_6B_umap_object_iso <- CreateSeuratObject(counts = G, project = "org_6B")
org_6C_umap_object_iso <- CreateSeuratObject(counts = H, project = "org_6C")


###
#merge the objects 
merged_seurat <- merge(org_1A_umap_object_iso, y = c(org_1B_umap_object_iso,
                                                org_3A_umap_object_iso, org_3B_umap_object_iso, org_3C_umap_object_iso, org_6A_umap_object_iso, org_6B_umap_object_iso, org_6C_umap_object_iso),
                       add.cell.ids = c("org1_A", "org1_B", "org3_A", "org3_B","org3_C", "org6_A", "org6_B", "org6_C"), project = 'cortex_org_iso')

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

## split sample column to makw a batch col
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Day', 'repliacte', 'Barcode'), 
                                    sep = '_')

## filter the data to iso and gene cells match
merged_seurat_isoform_filtered <- subset(merged_seurat, cells =obj@graphs[["RNA_nn"]]@Dimnames[[1]])

VlnPlot(merged_seurat_isoform_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_isoform_filtered <- NormalizeData(object = merged_seurat_isoform_filtered) # if using SCT dont run this
merged_seurat_isoform_filtered <- FindVariableFeatures(object = merged_seurat_isoform_filtered) # if using SCT dont run this
merged_seurat_isoform_filtered <- ScaleData(object = merged_seurat_isoform_filtered) # if using SCT dont run this
merged_seurat_isoform_filtered <- RunPCA(object = merged_seurat_isoform_filtered)
ElbowPlot(merged_seurat_isoform_filtered)
merged_seurat_isoform_filtered <- FindNeighbors(object = merged_seurat_isoform_filtered, dims = 1:10)
merged_seurat_isoform_filtered <- FindClusters(object = merged_seurat_isoform_filtered, resolution = 0.1)
merged_seurat_isoform_filtered <- FindClusters(object = merged_seurat_isoform_filtered, resolution = 0.3)
merged_seurat_isoform_filtered <- FindClusters(object = merged_seurat_isoform_filtered, resolution = 0.5)
merged_seurat_isoform_filtered <- FindClusters(object = merged_seurat_isoform_filtered, resolution = 0.7)
merged_seurat_isoform_filtered <- FindClusters(object = merged_seurat_isoform_filtered, resolution = 0.9)
merged_seurat_isoform_filtered <- RunUMAP(object = merged_seurat_isoform_filtered, dims = 1:10)

#plots
p1 <- DimPlot(merged_seurat_isoform_filtered, reduction = 'umap', group.by = 'orig.ident')
p2 <- DimPlot(merged_seurat_isoform_filtered, reduction = 'umap', group.by = 'Day')
p3 <- FeaturePlot(merged_seurat_isoform_filtered, reduction = 'umap', features = 'nCount_RNA')
p4 <- FeaturePlot(merged_seurat_isoform_filtered, reduction = "umap", features = 'nFeature_RNA')

pdf(file = "isoform_merged_objects.pdf", width = 12, height = 12) 
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off()

saveRDS(merged_seurat_isoform_filtered, file = "filt_seurat.merged.isofrom.rds")


######### ########

#### isofrom assay 
#Rejoin datasets after integration
merged_seurat_isoform_filtered <- JoinLayers(merged_seurat_isoform_filtered)
counts_table <- merged_seurat_isoform_filtered[["RNA"]]$counts

#obj <- integrated_harm_seurat

obj[["iso"]] <- CreateAssay5Object(counts = counts_table)

# Step 1: Normalize the new assay data
obj <- NormalizeData(obj, assay = "iso")
obj <- FindVariableFeatures(obj, assay = "iso")
obj <- ScaleData(obj, assay = "iso")

# Step 4: Perform PCA
obj <- RunPCA(obj, assay = "iso", reduction.name = "pca_iso")
# Step 5: Run UMAP
obj <- RunUMAP(obj, reduction = "pca_iso", dims = 1:10, assay = "iso", reduction.name = "umap_iso")

# Optionally, visualize the UMAP
DimPlot(obj, reduction = "umap_iso") | DimPlot(obj, reduction = "umap.harm")

#add in some metadata
#row.names(cell_groups) <- cell_groups$barcode_seq
#gene_counts$gene_id <- NULL

#pbmc <- AddMetaData(pbmc, metadata = cell_groups, col.name = 'groups')

FeaturePlot(obj, features = "VIM", reduction = "umap_iso") | FeaturePlot(obj, features = "VIM", reduction = "umap.harm")
DimPlot(obj, label = TRUE, reduction = "umap_iso") | DimPlot(obj, label = TRUE, reduction = "umap.harm")

DimPlot(obj, label = TRUE, reduction = "umap_iso", group.by = 'Day') 


saveRDS(obj, file = "no_bambu_gene_seurat_inter_harm_iso.rds")




