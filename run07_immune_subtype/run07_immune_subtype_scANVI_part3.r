library(Seurat)
library(SeuratWrappers)

setwd('/home/yq238/project/scRNA_seq/HTAN/07immune')

load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN2/06coarse_annotation/seurat_integrated_coarse_annotation.RData")
seurat_immune <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Immune'))


annotation <- read.csv('immune_cell_type.txt',sep='\t',header=FALSE)

seurat_immune@meta.data$immune_subtype = ""

for (i in 1:nrow(annotation)) {
    sample_id <- annotation[i, "V1"]
    subtype <- annotation[i, "V2"]
    if (sample_id %in% rownames(seurat_immune@meta.data)) {
        seurat_immune@meta.data[sample_id, "immune_subtype"] <- subtype
    }
}

seurat_immune <- FindVariableFeatures(seurat_immune)
seurat_immune <- ScaleData(seurat_immune)
seurat_immune <- RunPCA(seurat_immune)
seurat_immune <- IntegrateLayers(object = seurat_immune, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)

seurat_immune <- FindNeighbors(seurat_immune, reduction = "integrated.mnn", dims = 1:50)
seurat_immune <- FindClusters(seurat_immune, cluster.name = "mnn_clusters")
seurat_immune <- RunUMAP(seurat_immune, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")

pdf('dimplot_clusters_immune.pdf', width = 7 , height = 6)
DimPlot(
  seurat_immune,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()


pdf('immune_subtype_annotation_umap.pdf', width = 11 , height = 6)
DimPlot(seurat_immune, reduction = "umap.mnn", label = TRUE, repel = TRUE,group.by='immune_subtype')
dev.off()

save(seurat_immune, file = 'seurat_immune.RData')










