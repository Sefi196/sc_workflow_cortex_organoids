library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
library(presto)
library(harmony)

###
#merge the objects 
merged_seurat <- merge(org1A_umap_object, y = c(org1B_umap_object,
  org3A_umap_object, org3B_umap_object, org3C_umap_object, org6A_umap_object, org6B_umap_object, org6C_umap_object),
                       add.cell.ids = c("org1_A", "org1_B", "org3_A", "org3_B","org3_C", "org6_A", "org6_B", "org6_C"), project = 'cortex_org')

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

## split sample column to makw a batch col
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Day', 'repliacte', 'Barcode'), 
                                    sep = '_')
#check some QC metrics 
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Day")
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sctype_db")


table(merged_seurat$sctype_db)
table(merged_seurat$orig.ident)
table(merged_seurat$Day)

merged_seurat <- NormalizeData(object = merged_seurat) # if using SCT dont run this
merged_seurat <- FindVariableFeatures(object = merged_seurat) # if using SCT dont run this
merged_seurat <- ScaleData(object = merged_seurat) # if using SCT dont run this
merged_seurat <- RunPCA(object = merged_seurat)
ElbowPlot(merged_seurat)
merged_seurat <- FindNeighbors(object = merged_seurat, dims = 1:16)
merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.1)
merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.3)
merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.5)
merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.7)
merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.9)
merged_seurat <- RunUMAP(object = merged_seurat, dims = 1:30)


#plots
p1 <- DimPlot(merged_seurat, reduction = 'umap', group.by = 'orig.ident')
p2 <- DimPlot(merged_seurat, reduction = 'umap', group.by = 'sctype_db', label = T)
p3 <- FeaturePlot(merged_seurat, reduction = 'umap', features = 'nCount_RNA')
p4 <- FeaturePlot(merged_seurat, reduction = "umap", features = 'nFeature_RNA')
p5 <- DimPlot(merged_seurat, reduction = 'umap', group.by = 'Day')

pdf(file = "merged_objects.pdf", width = 12, height = 12) 
grid.arrange(p3, p4, p1, p5, p2, ncol = 2)
dev.off()

##save merged object
saveRDS(merged_seurat, file = "merged_seurat.rds")

merged_seurat <- JoinLayers(object = merged_seurat)
merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$orig.ident)


obj <- IntegrateLayers(
  object= merged_seurat,
  method= HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction="intergrated.harm",
  theta = 5,
  tau = 1,
  lambda = 1,
  npcs = 20,
  verbose=TRUE
)

obj <- FindNeighbors(obj, dims=1:30, reduction="intergrated.harm")
obj <- FindClusters(obj, resolution=0.7, cluster.name="harm_cluster")
obj <- RunUMAP(obj, reduction="intergrated.harm", dims=1:20, reduction.name = "umap.harm")

pdf(file = "integrated_harm.pdf", width = 12, height = 12) 
DimPlot(obj, reduction = "umap.harm", group.by = c("Day", "orig.ident", "harm_cluster", "sctype_db"), label = T)  
FeaturePlot(obj, reduction = "umap.harm", features = c('nFeature_RNA', "nCount_RNA"))
dev.off()

### save object 
saveRDS(obj, file = "integrated_harm_seurat.rds")


#### astrocytres 
FeaturePlot(obj, features = c("GFAP", "GAP43", "S100B", "ALDH1L1", "AQP4", "GLUL", "GJA1", "CD44", "SLC1A3"), reduction = "umap.harm")

#sctype GLUTA markers 
FeaturePlot(obj, features = c("SLC17A7","SLC17A6","GRIN1","GRIN2B","GLS","GLUL","GRIN2A"), reduction = "umap.harm")

#sctype DOPE makers 
FeaturePlot(obj, features = c("TH","SLC6A3","FOXA2","KCNJ6","NR4A2","LMX1B","DBH","SLC6A2","PPP1R1B"), reduction = "umap.harm")


FeaturePlot(obj, features = "GAD1", reduction = "umap.harm")


VlnPlot(obj, features = c("NR4A2"), split.by = "Day", group.by = "orig.ident")

## integrate with suerat 
obj <- IntegrateLayers(
  object= merged_seurat,
  method= CCAIntegration,
  orig.reduction = "pca",
  new.reduction="intergrated.CCA",
  verbose=TRUE
)

obj <- FindNeighbors(obj, dims=1:50, reduction="intergrated.CCA")
obj <- FindClusters(obj, resolution=0.7, cluster.name="cca_cluster")
obj <- RunUMAP(obj, reduction="intergrated.CCA", dims=1:20, reduction.name = "umap.cca")

pdf(file = "figures/integrated_cca.pdf", width = 12, height = 12) 
DimPlot(obj, reduction = "umap.cca", group.by = c("Day", "orig.ident", "cca_cluster", "sctype_db"), label = T)  
dev.off()


####
obj <- IntegrateLayers(
  object= merged_seurat,
  method= RPCAIntegration,
  orig.reduction = "pca",
  new.reduction="intergrated.RPCA",
  verbose=TRUE
)


obj <- FindNeighbors(obj, dims=1:50, reduction="intergrated.RPCA")
obj <- FindClusters(obj, resolution=0.7, cluster.name="rpca_cluster")
obj <- RunUMAP(obj, reduction="intergrated.RPCA", dims=1:10, reduction.name = "umap.RPCA")

pdf(file = "figures/integrated_rpca.pdf", width = 12, height = 12) 
DimPlot(obj, reduction = "umap.RPCA", group.by = c("Day", "orig.ident", "rpca_cluster", "sctype_db"), label = T)  
FeaturePlot(obj, reduction = "umap.RPCA", features = c('nFeature_RNA', "nCount_RNA"))
dev.off()

### look at the PCA
pdf(file = "figures/PCA.pdf", width = 6, height = 6) 
DimPlot(
  obj,
  reduction = "pca",
  group.by = "orig.ident",  # Specify the metadata column for individual samples
  label = TRUE,         # Optionally add labels to the clusters
  repel = TRUE          # Repel labels to avoid overlap
)
dev.off()


FeaturePlot(obj, features = c("HES1", "HOPX", "MOXD1"), reduction = "umap.cca")
FeaturePlot(obj, features = c("EGFR", "PDGFRA", "OLIG1", "OLIG2"), reduction = "umap.harm")

## small population of glial cells with molecular features characteristic of astroglia (GJA1, S100B, SPARC)."
FeaturePlot(obj, features = c("GJA1", "S100B", "SPARC"), reduction = "umap.cca")

##deep-layer projection neurons"
FeaturePlot(obj, features = c("FEZF1", "NEUROD2", "NEUROD6", "TBR1", "PDE1A"), reduction = "umap.harm")

#outer radial glia 
FeaturePlot(obj, features = c("SOX2", "HOPX", "PEA15", "LGALS3BP", "MOXD1"), reduction = "umap.harm")

#primarily apical radial glia (aRG; EMX1, SOX2, HES1) 
FeaturePlot(obj, features = c("EMX1", "SOX2", "HES1"), reduction = "umap.harm")

#intermediate progenitors (IP; EMX1, EOMES, INSM1), 
FeaturePlot(obj, features = c("EMX1", "EOMES", "INSM1"), reduction = "umap.harm")

#deep-layer corticofugal projection neurons (CFuPN; SOX5, LDB2, CRYM, TLE4)
FeaturePlot(obj, features = c("SOX5", "LDB2", "CRYM", "TLEA"), reduction = "umap.harm")

#Neurogenesis of callosal projection neurons (CPN; UNC5A, EPHA4, BHLHE22, PLXNA4, SATB2) 
FeaturePlot(obj, features = c("UNC5A", "EPHA4", "BHLHE22", "PLXNA4", "SATB2"), reduction = "umap.harm")



#Rejoin datasets after integration
seurat_integrated <- JoinLayers(obj)

#Idents(seurat_integrated) <- "sctype_db"
#### find Markers ###

#Find markers for every cluster compared to all remaining cells, report only the positive ones
cell.markers <- FindAllMarkers(seurat_integrated, assay = "RNA")

#Generate an expression heatmap plotting the top 5 upregulated markers for each cluster
cell.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n=3) %>%
  ungroup() -> top5

##### 
DoHeatmap(seurat_integrated, features=top5$gene)





######## set ideitites of cells ###### 

new_cluster_ids <- c('6' = 'Excitatory neuron',
                     '14' = 'Excitatory neuron', 
                     '5' = 'Excitatory neuron',
                     '16' = 'Excitatory neuron',
                     '3' = 'Excitatory neuron', 
                     '10' = 'ealy RG',
                     '11' = 'ealy RG',
                     '19' = 'ealy RG', 
                     '9' = 'ealy RG', 
                     '2' = 'Inhibitory neurons',
                     '15' = 'Neuronal Progenitor',
                     '1' = 'Neuronal Progenitor', 
                     '17' = 'Neuronal Progenitor',
                     '20' = 'Excitatory neuron',
                     '0' = 'Excitatory neuron',
                     '7' = 'oRG', 
                     '4' = 'Astrocytes',
                     '12' = 'Immature neurons', 
                     '13' = 'Immature neurons',
                     '8' = 'Immature neurons',
                     '18' = 'Immature neurons')

for_Mike_seurat_object <- RenameIdents(obj, new_cluster_ids)


pdf(file = "figure_for_mike_ivestigator_grant.pdf", width = 6, height = 6) 
DimPlot(for_Mike_seurat_object, reduction = "umap.harm", label = T)  
dev.off()


FeaturePlot(obj, features = c("MAPT", "GRIA1", "BHLHE22", "PLXNA4", "SATB2"), reduction = "umap.harm")


