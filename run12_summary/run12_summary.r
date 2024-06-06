library(Seurat)

setwd('/home/yq238/project/scRNA_seq/HTAN/12summary')

load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN/06coarse_annotation/seurat_integrated_coarse_annotation.RData")
load('/home/yq238/project/scRNA_seq/HTAN/07immune/seurat_immune.RData')
load('/home/yq238/project/scRNA_seq/HTAN/09SCLC_NSCLC/seurat_NSCLC_join.RData')
load('/home/yq238/project/scRNA_seq/HTAN/10SCLC_subtype/seurat_SCLC_join_annotation.RData')
load('/home/yq238/project/scRNA_seq/HTAN/11normal_subtype/seurat_normal.RData')

seurat_integrated@meta.data$final_annotation = ""


for (i in 1:nrow(seurat_immune@meta.data)) {
    sample_id <- rownames(seurat_immune@meta.data)[i]
    subtype <- seurat_immune@meta.data$immune_subtype[i]
    seurat_integrated@meta.data[sample_id, "final_annotation"] <- subtype
}

for (i in 1:nrow(seurat_NSCLC_join@meta.data)) {
    sample_id <- rownames(seurat_NSCLC_join@meta.data)[i]
    subtype <- 'NSCLC'
    seurat_integrated@meta.data[sample_id, "final_annotation"] <- subtype
}

for (i in 1:nrow(seurat_SCLC_join@meta.data)) {
    sample_id <- rownames(seurat_SCLC_join@meta.data)[i]
    subtype <- seurat_SCLC_join@meta.data$SCLC_subtype[i]
    seurat_integrated@meta.data[sample_id, "final_annotation"] <- subtype
}

for (i in 1:nrow(seurat_normal@meta.data)) {
    sample_id <- rownames(seurat_normal@meta.data)[i]
    subtype <- seurat_normal@meta.data$epithelial_normal_subtype[i]
    seurat_integrated@meta.data[sample_id, "final_annotation"] <- subtype
}

Mesenchymal_cells <- rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$coarse_annotation == "Mesenchymal",])
seurat_integrated@meta.data$final_annotation[rownames(seurat_integrated@meta.data) %in% Mesenchymal_cells] <- 'Mesenchymal'



pdf('final_annotation_plot.pdf', width = 14 , height = 10)
DimPlot(seurat_integrated, reduction = "umap.mnn", label = TRUE, repel = TRUE, group.by='final_annotation')
dev.off()



save(seurat_integrated, file='seurat_integrated_final_annotation.RData')




