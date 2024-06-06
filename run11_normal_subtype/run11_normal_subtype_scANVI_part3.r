library(Seurat)
library(SeuratWrappers)

setwd('/home/yq238/project/scRNA_seq/HTAN/11normal_subtype')

load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN/06coarse_annotation/seurat_integrated_coarse_annotation.RData")

seurat_epithelial <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Epithelial'))

tumor_cells <- read.table("../08SCLC_NSCLC/cancer_cells.txt")
tumor_cells <- tumor_cells$V1

seurat_epithelial@meta.data$tumorstatus <- ifelse(rownames(seurat_epithelial@meta.data) %in% tumor_cells, "cancer",
 "normal")

seurat_normal <- subset(x = seurat_epithelial, subset = (tumorstatus == 'normal'))


annotation <- read.csv('normal_cell_type.txt',sep='\t',header=FALSE)

seurat_normal@meta.data$epithelial_normal_subtype = ""


for (i in 1:nrow(annotation)) {
    sample_id <- annotation[i, "V1"]
    subtype <- annotation[i, "V2"]
    if (sample_id %in% rownames(seurat_normal@meta.data)) {
        seurat_normal@meta.data[sample_id, "epithelial_normal_subtype"] <- subtype
    }
}

seurat_normal <- FindVariableFeatures(seurat_normal)
seurat_normal <- ScaleData(seurat_normal)
seurat_normal <- RunPCA(seurat_normal)
seurat_normal <- IntegrateLayers(object = seurat_normal, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)

seurat_normal <- FindNeighbors(seurat_normal, reduction = "integrated.mnn", dims = 1:50)
seurat_normal <- FindClusters(seurat_normal, cluster.name = "mnn_clusters")
seurat_normal <- RunUMAP(seurat_normal, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")


pdf('normal_subtype_annotation_umap.pdf', width = 11 , height = 6)
DimPlot(seurat_normal, reduction = "umap.mnn", label = TRUE, repel = TRUE,group.by='epithelial_normal_subtype')
dev.off()

save(seurat_normal, file = 'seurat_normal.RData')










