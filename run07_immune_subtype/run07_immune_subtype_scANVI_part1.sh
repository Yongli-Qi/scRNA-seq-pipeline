library(Seurat)

setwd('/home/yq238/project/scRNA_seq/HTAN/07immune')

load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN/06coarse_annotation/seurat_integrated_coarse_annotation.RData")
seurat_immune <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Immune'))


seurat_immune_join <- JoinLayers(seurat_immune)

sce_obj <- as.SingleCellExperiment(seurat_immune_join, assay = c("RNA"))
library(zellkonverter)
writeH5AD(sce_obj, "sce_obj.h5ad", X_name = 'counts')
