library(Seurat)

setwd('/home/yq238/project/scRNA_seq/HTAN/11normal_subtype')

load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN/06coarse_annotation/seurat_integrated_coarse_annotation.RData")

seurat_epithelial <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Epithelial'))

tumor_cells <- read.table("../09SCLC_NSCLC/cancer_cells.txt")
tumor_cells <- tumor_cells$V1

seurat_epithelial@meta.data$tumorstatus <- ifelse(rownames(seurat_epithelial@meta.data) %in% tumor_cells, "cancer",
 "normal")

seurat_normal <- subset(x = seurat_epithelial, subset = (tumorstatus == 'normal'))


seurat_normal_join <- JoinLayers(seurat_normal)

sce_obj <- as.SingleCellExperiment(seurat_normal_join, assay = c("RNA"))
library(zellkonverter)
writeH5AD(sce_obj, "sce_obj.h5ad", X_name = 'counts')

