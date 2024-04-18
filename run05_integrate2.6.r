library(Seurat)
library(SeuratWrappers)

setwd("/home/yq238/project/scRNA_seq/HTAN2/05seurat_obj_integrated")
load("../04additional_check/seurat_combined_PCA.RData")

options(future.globals.maxSize = 10000 * 1024^2)

seurat_integrated <- IntegrateLayers(object = seurat_combined, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)
save(seurat_integrated, file = 'seurat_integrated_fastmnn.RData')

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "integrated.mnn", dims = 1:50)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 2.6, cluster.name = "mnn_clusters")

seurat_integrated <- RunUMAP(seurat_integrated, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
pdf('dimplot_integrated_mnn.pdf', width = 7 , height = 6)
DimPlot(
  seurat_integrated ,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()


# output the cell names
write.table(rownames(seurat_integrated@meta.data),'cells_mine.txt',sep = '\t',quote = F, col.names = F, row.names = F)

