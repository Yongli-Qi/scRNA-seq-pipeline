library(Seurat)


setwd('/home/yq238/project/scRNA_seq/HTAN/09SCLC_NSCLC')
load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN/06coarse_annotation/seurat_integrated_coarse_annotation.RData")

seurat_epithelial <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Epithelial'))

tumor_cells <- read.table("cancer_cells.txt")
tumor_cells <- tumor_cells$V1

seurat_epithelial@meta.data$tumorstatus <- ifelse(rownames(seurat_epithelial@meta.data) %in% tumor_cells, "cancer",
 "normal")

pdf('dimplot_clusters_tumor_normal.pdf', width = 7 , height = 6)
DimPlot(
  seurat_epithelial,
  reduction = "umap.mnn",
  group.by = c("tumorstatus"),
  combine = FALSE, label.size = 2
)
dev.off()



seurat_epithelial_join <- JoinLayers(seurat_epithelial)

seurat_tumor_join <- subset(x = seurat_epithelial_join, subset = (tumorstatus == 'cancer'))

seurat_tumor_join@meta.data$sample = ""
seurat_tumor_join@meta.data$sample <- ifelse(grepl("^HTA8_1", seurat_tumor_join@meta.data$orig.ident), "NSCLC", 
                           ifelse(grepl("^HTA8_2", seurat_tumor_join@meta.data$orig.ident), "SCLC", seurat_tumor_join@meta.data$sample))

seurat_SCLC_join <- subset(x = seurat_tumor_join, subset = (sample == 'SCLC'))
seurat_NSCLC_join <- subset(x = seurat_tumor_join, subset = (sample == 'NSCLC'))


pdf('dimplot_clusters_tumor_SCLC.pdf', width = 7 , height = 6)
DimPlot(
  seurat_SCLC_join,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()

pdf('dimplot_clusters_tumor_NSCLC.pdf', width = 7 , height = 6)
DimPlot(
  seurat_NSCLC_join,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()



pdf('tumor_subtype_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_tumor_join, reduction = "umap.mnn", label = TRUE, repel = TRUE,group.by='sample')
dev.off()



pdf('gene_expression_violin_CHGA.pdf', width = 5 , height = 6)
VlnPlot(seurat_tumor_join, features = c("CHGA"),group.by = 'sample',pt.size = 0)
dev.off()
pdf('gene_expression_violin_CHGB.pdf', width = 5 , height = 6)
VlnPlot(seurat_tumor_join, features = c("CHGB"),group.by = 'sample',pt.size = 0)
dev.off()
pdf('gene_expression_violin_SYP.pdf', width = 5 , height = 6)
VlnPlot(seurat_tumor_join, features = c("SYP"),group.by = 'sample',pt.size = 0)
dev.off()
pdf('gene_expression_violin_ASCL1.pdf', width = 5 , height = 6)
VlnPlot(seurat_tumor_join, features = c("ASCL1"),group.by = 'sample',pt.size = 0)
dev.off()
pdf('gene_expression_violin_INSM1.pdf', width = 5 , height = 6)
VlnPlot(seurat_tumor_join, features = c("INSM1"),group.by = 'sample',pt.size = 0)
dev.off()
pdf('gene_expression_violin_BEX1.pdf', width = 5 , height = 6)
VlnPlot(seurat_tumor_join, features = c("BEX1"),group.by = 'sample',pt.size = 0)
dev.off()


save(seurat_SCLC_join,file='seurat_SCLC_join.RData')
save(seurat_NSCLC_join,file='seurat_NSCLC_join.RData')





