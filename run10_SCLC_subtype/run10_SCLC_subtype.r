library(Seurat)


setwd("/home/yq238/project/scRNA_seq/HTAN/10SCLC_subtype")
load("/home/yq238/project/scRNA_seq/HTAN/09SCLC_NSCLC/seurat_SCLC_join.RData")


seurat_SCLC_join <- FindVariableFeatures(seurat_SCLC_join)
seurat_SCLC_join <- ScaleData(seurat_SCLC_join)
seurat_SCLC_join <- RunPCA(seurat_SCLC_join)

save(seurat_SCLC_join, file="seurat_SCLC_join_rePCA.RData")


total_variance <- seurat_SCLC_join@reductions$pca@misc$total.variance
eigValues = (seurat_SCLC_join[["pca"]]@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)

seurat_SCLC_join <- FindNeighbors(seurat_SCLC_join, reduction = "pca", dims = 1:50)
seurat_SCLC_join <- FindClusters(seurat_SCLC_join, resolution = 0.8, cluster.name = "clusters")
seurat_SCLC_join <- RunUMAP(seurat_SCLC_join, reduction = "pca", dims = 1:50, reduction.name = "umap")


pdf('dimplot_clusters_sclc_0.8.pdf', width = 7 , height = 6)
DimPlot(
  seurat_SCLC_join,
  reduction = "umap",
  group.by = c("clusters"),
  combine = FALSE, label.size = 2
)
dev.off()





########################################

pdf('gene_expression_violin_SCLC_ASCL1.pdf', width = 14, height = 6)
VlnPlot(seurat_SCLC_join, features = c("ASCL1"),pt.size = 0)
dev.off()


pdf('gene_expression_violin_SCLC_NEUROD1.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("NEUROD1"),pt.size = 0)
dev.off()

pdf('gene_expression_violin_SCLC_YAP1.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("YAP1"),pt.size = 0)
dev.off()

pdf('gene_expression_violin_SCLC_POU2F3.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("POU2F3"),pt.size = 0)
dev.off()


pdf('gene_expression_violin_SCLC_CHGA.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("CHGA"),pt.size = 0)
dev.off()

pdf('gene_expression_violin_SCLC_CHGB.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("CHGB"),pt.size = 0)
dev.off()

pdf('gene_expression_violin_SCLC_SYP.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("SYP"),pt.size = 0)
dev.off()

pdf('gene_expression_violin_SCLC_BEX1.pdf', width = 14, height = 6)
VlnPlot(seurat_SCLC_join, features = c("BEX1"),pt.size = 0)
dev.off()

pdf('gene_expression_violin_SCLC_KDM1A.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("KDM1A"),pt.size = 0)
dev.off()

pdf('gene_expression_violin_SCLC_INSM1.pdf', width = 14 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("INSM1"),pt.size = 0)
dev.off()


######################################################################

seurat_SCLC_join@meta.data$SCLC_subtype <- ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(0,5,6,7,9,10,11,13,14,15,17,19,20,22,24,25,26,27,29,30,32,33), "SCLC-A",
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(2,3,4,16,31,34), "SCLC-N", 
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(8,18), "SCLC-P", 
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(1,21,23), "SCLC-NA-NE",
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(12,28), "SCLC-NA-nonNE",
"NA")))))

pdf('SCLC_subtype_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_SCLC_join, reduction = "umap", label = TRUE, repel = TRUE, group.by='SCLC_subtype')
dev.off()

table(seurat_SCLC_join@meta.data$SCLC_subtype)

#SCLC-A        SCLC-N    SCLC-NA-NE SCLC-NA-nonNE        SCLC-P 
#        29040         12848          5701          1951          2951 

#######output
SCLCA_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-A",])
write(SCLCA_cells, file = "SCLC-A_cells.txt")

SCLCN_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-N",])
write(SCLCN_cells, file = "SCLC-N_cells.txt")

SCLCP_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-P",])
write(SCLCP_cells, file = "SCLC-P_cells.txt")

SCLC_NA_NE_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-NA-NE",])
write(SCLC_NA_NE_cells, file = "SCLC-NA-NE_cells.txt")

SCLC_NA_nonNE_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-NA-nonNE",])
write(SCLC_NA_nonNE_cells, file = "SCLC_NA_nonNE_cells.txt")

save(seurat_SCLC_join,file='seurat_SCLC_join_annotation.RData')
