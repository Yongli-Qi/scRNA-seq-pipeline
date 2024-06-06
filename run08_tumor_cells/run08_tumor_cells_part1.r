library(Seurat)


setwd('/home/yq238/project/scRNA_seq/HTAN/08tumor_cells')

load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN/06coarse_annotation/seurat_integrated_coarse_annotation.RData")

seurat_epithelial <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Epithelial'))

seurat_epithelial_join <- JoinLayers(seurat_epithelial)
count_matrix <- GetAssayData(seurat_epithelial_join, layer = "counts")
#save(count_matrix, file="count_matrix_epithelial_sparse.RData")

count_matrix <- as.matrix(count_matrix)
save(count_matrix, file="count_matrix_epithelial.RData")

