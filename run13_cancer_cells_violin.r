library(Seurat)

setwd("/home/yq238/project/scRNA_seq/HTAN2/13cancer_cells_violin")
load("../13cancer_cells_semi_supervised/seurat_SCLC_join.RData")

seurat_SCLC_join <- FindVariableFeatures(seurat_SCLC_join)
seurat_SCLC_join <- ScaleData(seurat_SCLC_join)
seurat_SCLC_join <- RunPCA(seurat_SCLC_join)

save(seurat_SCLC_join, file="seurat_SCLC_join_rePCA.RData")

seurat_SCLC_join <- FindNeighbors(seurat_SCLC_join, reduction = "pca", dims = 1:50)
seurat_SCLC_join <- FindClusters(seurat_SCLC_join, resolution = 0.4, cluster.name = "clusters")
seurat_SCLC_join <- RunUMAP(seurat_SCLC_join, reduction = "pca", dims = 1:50, reduction.name = "umap")

pdf('dimplot_clusters_sclc_0.4.pdf', width = 7 , height = 6)
DimPlot(
  seurat_SCLC_join,
  reduction = "umap",
  group.by = c("clusters"),
  combine = FALSE, label.size = 2
)
dev.off()


########################################

pdf('gene_expression_violin_ASCL1.pdf', width = 10, height = 6)
VlnPlot(seurat_SCLC_join, features = c("ASCL1"),split.by = 'seurat_clusters')
dev.off()


pdf('gene_expression_violin_NEUROD1.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("NEUROD1"),split.by = 'seurat_clusters')
dev.off()

pdf('gene_expression_violin_YAP1.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("YAP1"),split.by = 'seurat_clusters')
dev.off()

pdf('gene_expression_violin_POU2F3.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("POU2F3"),split.by = 'seurat_clusters')
dev.off()


pdf('gene_expression_violin_CHGA.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("CHGA"),split.by = 'seurat_clusters')
dev.off()

pdf('gene_expression_violin_CHGB.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("CHGB"),split.by = 'seurat_clusters')
dev.off()

pdf('gene_expression_violin_SYP.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("SYP"),split.by = 'seurat_clusters')
dev.off()

pdf('gene_expression_violin_BEX1.pdf', width = 10, height = 6)
VlnPlot(seurat_SCLC_join, features = c("BEX1"),split.by = 'seurat_clusters')
dev.off()

pdf('gene_expression_violin_KDM1A.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("KDM1A"),split.by = 'seurat_clusters')
dev.off()

pdf('gene_expression_violin_INSM1.pdf', width = 10 , height = 6)
VlnPlot(seurat_SCLC_join, features = c("INSM1"),split.by = 'seurat_clusters')
dev.off()

######################################################################

seurat_SCLC_join@meta.data$SCLC_subtype <- ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(0,5,6,7,8,9,10,11,13,14,16,20,21,22), "SCLC-A",
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(2,3,4,15), "SCLC-N", 
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(12), "SCLC-P", 
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(18,19), "SCLC-Y",
"NA"))))

pdf('SCLC_subtype_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_SCLC_join, reduction = "umap", label = TRUE, repel = TRUE, group.by='SCLC_subtype')
dev.off()

table(seurat_SCLC_join@meta.data$SCLC_subtype)

 #   NA SCLC-A SCLC-N SCLC-P SCLC-Y 
 # 5066  28409  12035   1335   1183


#######output
SCLCA_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-A",])
write(SCLCA_cells, file = "SCLC-A_cells.txt")

SCLCN_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-N",])
write(SCLCN_cells, file = "SCLC-N_cells.txt")

SCLCP_cells <- rownames(seurat_SCLC_join@meta.data[seurat_SCLC_join@meta.data$SCLC_subtype == "SCLC-P",])
write(SCLCP_cells, file = "SCLC-P_cells.txt")






