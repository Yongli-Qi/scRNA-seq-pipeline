library(Seurat)
.libPaths(c("/vast/palmer/apps/avx2/software/R/4.3.0-foss-2022b/lib64/R/library", .libPaths()))
library(SeuratWrappers)

setwd('/home/yq238/project/scRNA_seq/HTAN2/12extract_NSCLC_SCLC')
load("../10_epithelial_tumors_raw_counts/seurat_epithelial_tumorstatus_annotation.RData")


seurat_tumor <- subset(x = seurat_epithelial, subset = (tumorstatus == 'cancer'))


############# add annotation
seurat_tumor@meta.data$cancer_type <- ifelse(seurat_tumor@meta.data$ne_subtype == "non-neuroendocrine", "NSCLC", "SCLC")


############## extract the cells
SCLC_cells <- rownames(seurat_tumor@meta.data[seurat_tumor@meta.data$cancer_type == "SCLC",])
write(SCLC_cells, file = "SCLC cell.txt")

NSCLC_cells <- rownames(seurat_tumor@meta.data[seurat_tumor@meta.data$cancer_type == "NSCLC",])
write(NSCLC_cells, file = "NSCLC cell.txt")

############plot
pdf('cancer_cells_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_tumor, reduction = "umap.mnn", label = TRUE, repel = TRUE, group.by='cancer_type')
dev.off()

##########save

seurat_SCLC <- subset(x = seurat_tumor, subset = (cancer_type == 'SCLC'))
#seurat_NSCLC <- subset(x = seurat_tumor, subset = (cancer_type == 'NSCLC'))

save(seurat_SCLC, file="seurat_SCLC.RData")
#save(seurat_NSCLC, file="seurat_NSCLC.RData")








