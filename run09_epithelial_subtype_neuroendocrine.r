library(Seurat)
library(SeuratWrappers)

setwd('/home/yq238/project/scRNA_seq/HTAN2/09_epithelial_subtype_revised')
load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN2/06coarse_annotation/seurat_integrated_coarse_annotation.RData")

seurat_epithelial <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Epithelial'))

seurat_epithelial <- FindVariableFeatures(seurat_epithelial)
seurat_epithelial <- ScaleData(seurat_epithelial)
seurat_epithelial <- RunPCA(seurat_epithelial)
seurat_epithelial <- IntegrateLayers(object = seurat_epithelial, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)

save(seurat_epithelial, file="seurat_epithelial_integrated.RData")

seurat_epithelial <- FindNeighbors(seurat_epithelial, reduction = "integrated.mnn", dims = 1:50)
seurat_epithelial <- FindClusters(seurat_epithelial, resolution = 2, cluster.name = "mnn_clusters")
seurat_epithelial <- RunUMAP(seurat_epithelial, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")


pdf('dimplot_clusters_epithelial_mnn_2.pdf', width = 7 , height = 6)
DimPlot(
  seurat_epithelial,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()


total_variance <- seurat_epithelial@reductions$pca@misc$total.variance
eigValues = (seurat_epithelial[["pca"]]@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)


seurat_epithelial_join <- JoinLayers(seurat_epithelial)

pdf('gene_expression_umap_3markers_2.pdf', width = 7 , height = 6)
FeaturePlot(seurat_epithelial_join, 
             reduction = "umap.mnn", 
             features = c("CHGA","CHGB","SYP"), 
             order = TRUE,
             min.cutoff = 'q10', 
             label = TRUE, raster=FALSE)
dev.off()


pdf('gene_expression_violin_CHGA.pdf', width = 14 , height = 6)
VlnPlot(seurat_epithelial_join, features = c("CHGA"),split.by = 'seurat_clusters')
dev.off()
pdf('gene_expression_violin_CHGB.pdf', width = 14 , height = 6)
VlnPlot(seurat_epithelial_join, features = c("CHGB"),split.by = 'seurat_clusters')
dev.off()
pdf('gene_expression_violin_SYP.pdf', width = 14 , height = 6)
VlnPlot(seurat_epithelial_join, features = c("SYP"),split.by = 'seurat_clusters')
dev.off()
pdf('gene_expression_violin_ASCL1.pdf', width = 14 , height = 6)
VlnPlot(seurat_epithelial_join, features = c("ASCL1"),split.by = 'seurat_clusters')
dev.off()
pdf('gene_expression_violin_INSM1.pdf', width = 14 , height = 6)
VlnPlot(seurat_epithelial_join, features = c("INSM1"),split.by = 'seurat_clusters')
dev.off()
pdf('gene_expression_violin_BEX1.pdf', width = 14 , height = 6)
VlnPlot(seurat_epithelial_join, features = c("BEX1"),split.by = 'seurat_clusters')
dev.off()


seurat_epithelial@meta.data$ne_subtype <- ifelse(seurat_epithelial@meta.data$seurat_cluster %in% c(3,4,5,6,7,9,10,11,12,13,14,15,16,17,20,22,24,25,26,27,29,30,31,32,33,35,36,38,41), "neuroendocrine","non-neuroendocrine")

pdf('ne_cells_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_epithelial, reduction = "umap.mnn", label = TRUE, repel = TRUE, group.by='ne_subtype')
dev.off()


############################################################
#        extract the names of cells in each type
############################################################

ne_cells <- rownames(seurat_epithelial@meta.data[seurat_epithelial@meta.data$ne_subtype == "neuroendocrine",])
non_ne_cells <- rownames(seurat_epithelial@meta.data[seurat_epithelial@meta.data$ne_subtype == "non-neuroendocrine",])

write(ne_cells, file = "ne_cell.txt")
write(non_ne_cells, file = "non-ne_cell.txt")


###ne cells
ref <- readRDS('ref_epi.rds')
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Neuroendocrine',]),'Neuroendocrine_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='SCLC-A',]),'SCLC-A_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='SCLC-N',]),'SCLC-N_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='SCLC-P',]),'SCLC-P_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)


####non_ne
selected_cells <- rownames(ref@meta.data[!ref@meta.data$author_cell_type %in% c("SCLC-P", "SCLC-A", "SCLC-N", "Neuroendocrine"),])
write.table(selected_cells, 'non_Neuroendocrine_cells_ref.txt', sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)


save(seurat_epithelial, file="seurat_epithelial_ne_annotation.RData")
