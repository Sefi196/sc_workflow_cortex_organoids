library(dplyr)
library(Seurat)
library(UpSetR)

##### This is unrelrated but i tried to do Upset plot for gene and isofrms that are uniqe to each cluster. 
### seems to work ok but its not so interesing as nearly all gene are found in all cell lines 
obj <- seurat_inter_harm_iso
JoinLayers(obj) -> obj
# Extract cluster information
cluster_info <- obj$seurat_clusters

# Get expression matrix and calculate gene cluster summary directly
iso_cluster_summary <- GetAssayData(obj, assay = "iso", layer = "counts") %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  tidyr::gather(key = "cell", value = "expression", -gene) %>%
  mutate(cluster = cluster_info[cell]) %>%
  group_by(gene, cluster) %>%
  summarize(total_expression = sum(expression, na.rm = TRUE), .groups = 'drop') %>%
  filter(total_expression > 0) %>%
  group_by(cluster) %>%
  summarize(genes = list(gene), .groups = 'drop')

cluster_gene_list <- setNames(gene_cluster_summary$genes, gene_cluster_summary$cluster)
cluster_iso_list <- setNames(iso_cluster_summary$genes, iso_cluster_summary$cluster)

#to plot them all 
upset(fromList(cluster_gene_list), nsets = length(cluster_gene_list), nintersects = NA)


##extact the sets you want
# Extract only the sets for clusters '0' and '4'
specific_clusters <- c('0', '4', "3")
set_order <- c('0', '4', '3')

# Filter the list to keep only the specified clusters
filtered_gene_list <- cluster_gene_list[specific_clusters]
filtered_iso_list <- cluster_iso_list[specific_clusters]

G <- upset(fromList(filtered_gene_list), sets = set_order, 
           #order.by = "freq",   # Or use "degree" or "matrix" depending on your preference
      keep.order = TRUE)
I <- upset(fromList(filtered_iso_list), sets = set_order, 
           #order.by = "freq",   # Or use "degree" or "matrix" depending on your preference
      keep.order = TRUE) 

pdf("UpsetPlots_genes_isoforms_in each cluster", width = 12, height = 6)
G + tit("Unique Gene in each cluster",x = 0.65, y=0.95, gp=gpar(fontsize=20))
I + grid.text("Unique Isoforms in each cluster",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

saveRDS()

DimPlot(obj, reduction = "umap.harm", label = T)


