library(Matrix)
library(dplyr)
library(rtracklayer)
library(Seurat)

# Function to find gene ID and transcripts for a given gene symbol
findGeneTranscripts <- function(geneSymbol, mappingFilePath, countMatrix) {
  # Read the mapping file to get gene ID for the provided gene symbol
  mapping <- read.csv(mappingFilePath, header = TRUE, stringsAsFactors = FALSE)
  geneID <- mapping %>% 
    filter(GeneSymbol == geneSymbol) %>% 
    pull(gene_id)
  
  # Read the count matrix from the provided file path
  countMatrix <- read.csv(countMatrixPath, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract the transcript IDs for the given gene ID
  transcripts <- unique(countMatrix[countMatrix[, 2] == geneID, 1])
  
  # Return the gene ID and transcript IDs
  return(list(GeneID = geneID, Transcripts = transcripts))
}

###### Running the fucntion 

# Mapping file path (GeneSymbol to GeneID)
mappingFilePath <- "/data/projects/punim1441/FLAMES_202311/resources/mapping_file_find_gene_transcript.csv"


# Mapping file path (GeneSymbol to GeneID)
###file for 20_min_sup_cnt
#countMatrixPath <- "/data/gpfs/projects/punim1441/Project_Kolf_SCLR/rebase/flames_all_together_supcnt_20/gene_isoform_seurat_objects/geneid_ENSID_raw.all.transcript.counts.csv"
#countMatrixPath <- "/data/gpfs/projects/punim1441/Project_Kolf_SCLR/flames_rebase_all_together/bigman/min_sup_cnt_100/min_sup_cnt100_geneid_ENSID_raw.all.transcript.counts.csv"
#countMatrixPath <- "/data/projects/punim1441/yairp/FLAMES_202311/resources/NDRD_gene_transcript_isd.csv"
countMatrixPath <- "gene_transcript_isd.csv"
#

# Find gene ID and transcript IDs for the given gene symbol

geneSymbol <- "MAPT"

result <- findGeneTranscripts(geneSymbol, mappingFilePath, countMatrixPath)
#if you want to use them for ploting with feature plot as it stands we need to change _ to -
result$Transcripts <- gsub("_", "-", result$Transcripts) %>% paste0("-", geneSymbol)

###can print out the 
# Print the gene ID and transcript IDs
cat("Gene ID for", geneSymbol, "is:", result$GeneID, "\n")
cat("Transcript IDs for", geneSymbol, "are:", result$Transcripts, "\n")

FeaturePlot(obj, features = geneSymbol, reduction = "umap.harm")
FeaturePlot(obj, features = result$Transcripts, reduction = "umap.harm")


DimPlot(topup_merged_integrated.scanorama, group.by =  "Day") | DimPlot(Isoform_topup_merged_integrated.scanorama, group.by =  "Day") 
#FeaturePlot(minsupcnt_20_seurat.integrated.isofrom, features = result$Transcripts)
dev.off()

##### PKM ########
result$Transcripts

pdf("Isoforms_CLU_wide.pdf", width = 6, height = 3)
FeaturePlot(obj, features = c("ENST00000405140.7"), reduction = "umap.harm")
dev.off()


pdf("Isoforms_GPM6B_UMAP.pdf", width = 15, height = 3)
FeaturePlot(Isoform_topup_merged_integrated.scanorama,
            features = c("ENST00000454189.6","BambuTx4994", "ENST00000316715.9", "ENST00000355135.6"),
            ncol = 4, pt.size = 0.2,
            keep.scale = "all", min.cutoff = 'q10')
dev.off()



pdf("Gene_PKM_UMAP.pdf", width = 5, height = 3)
FeaturePlot(topup_merged_integrated.scanorama,
            features = c("PKM"), pt.size = 0.2)
dev.off()


pdf("Gene_isoforms_byDay_UMAP.pdf", width = 6, height = 6)
DimPlot(topup_merged_integrated.scanorama, group.by =  "Day") 
DimPlot(Isoform_topup_merged_integrated.scanorama, group.by =  "Day") 
dev.off()



VlnPlot(object = Isoform_topup_merged_integrated.scanorama,
        #idents  = c("Mature neurons", "Neural Progenitor cells", "Immature neurons"),
        features = c("ENST00000335181.10", "ENST00000389093.7", "BambuTx3706"),  same.y.lims = TRUE)


topup_merged_integrated.scanorama$Day <- factor(topup_merged_integrated.scanorama$Day, levels = c("STC", "Day25", "Day55", "Day80"))
Isoform_topup_merged_integrated.scanorama$Day <- factor(Isoform_topup_merged_integrated.scanorama$Day, levels = c("STC", "Day25", "Day55", "Day80"))


#### set ident to plot subsets of cell types ####
Isoform_topup_merged_integrated.scanorama <- SetIdent(Isoform_topup_merged_integrated.scanorama, value = "prelim_cell_types")
Isoform_topup_merged_integrated.scanorama@active.ident <- factor(x = Isoform_topup_merged_integrated.scanorama@active.ident, levels = c("Stem Cell", "Radial Glia Early", "GABAergic interneuron Cluster 1"))


pdf("Isoforms_PKM_by_celltype.pdf", width = 8, height = 6)
# Now create the violin plot
VlnPlot(
  object = Isoform_topup_merged_integrated.scanorama,
  features = c("ENST00000335181.10", "ENST00000389093.7", "BambuTx3706"),
  #group.by = "prelim_cell_types",
  idents  = c("Stem Cell", "Radial Glia Early", "GABAergic interneuron Cluster 1"),
  same.y.lims = TRUE,
  assay = "RNA",
  slot = 'data'
)
dev.off()

DotPlot(object = Isoform_topup_merged_integrated.scanorama, features = c("BambuTx3706"), group.by = "orig.ident")

pdf("Gene_PKM_by_day.pdf", width = 5, height = 3)
VlnPlot(object =  topup_merged_integrated.scanorama,
        features = c("PKM"),
        group.by = "Day", same.y.lims = TRUE, assay="RNA", slot = 'data')
dev.off()



DE_isofrms <-read.csv("isoform_with_cell_types_markers.all_fc_1.clusters.csv")
DoHeatmap(object = Isoform_topup_merged_integrated.scanorama, features = DE_isofrms$gene, group.by = "ident")
