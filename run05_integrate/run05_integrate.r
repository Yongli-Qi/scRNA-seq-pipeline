library(Seurat)
library(SeuratWrappers)

setwd("/home/yq238/project/scRNA_seq/HTAN/05seurat_obj_integrated")
load("../04additional_check/seurat_combined_PCA.RData")

options(future.globals.maxSize = 10000 * 1024^2)

seurat_integrated <- IntegrateLayers(object = seurat_combined, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)
save(seurat_integrated, file = 'seurat_integrated_fastmnn.RData')



# output the cell names
write.table(rownames(seurat_integrated@meta.data),'cells_mine.txt',sep = '\t',quote = F, col.names = F, row.names = F)
