library(Seurat)
library(SeuratWrappers)

setwd('/home/yq238/project/scRNA_seq/HTAN2/07immune_subtype_revised')
load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN2/06coarse_annotation/seurat_integrated_coarse_annotation.RData")

seurat_immune <- subset(x = seurat_integrated, subset = (coarse_annotation == 'Immune'))

seurat_immune <- FindVariableFeatures(seurat_immune)
seurat_immune <- ScaleData(seurat_immune)
seurat_immune <- RunPCA(seurat_immune)
seurat_immune <- IntegrateLayers(object = seurat_immune, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)

save(seurat_immune, file="seurat_immune_integrated.RData")

seurat_immune <- FindNeighbors(seurat_immune, reduction = "integrated.mnn", dims = 1:50)
seurat_immune <- FindClusters(seurat_immune, resolution = 0.4, cluster.name = "mnn_clusters")
seurat_immune <- RunUMAP(seurat_immune, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")

pdf('dimplot_clusters_immune_mnn.pdf', width = 7 , height = 6)
DimPlot(
  seurat_immune,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()


total_variance <- seurat_immune@reductions$pca@misc$total.variance
eigValues = (seurat_immune[["pca"]]@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)


lines <- readLines('Data S2_combined.gmt')

cellTypeToMarkers <- list()
for(line in lines) {
  elements <- strsplit(line, "\t")[[1]]
  cellType <- elements[1]
  markers <- elements[-(1:2)]
  markers <- markers[markers != ""]
  cellTypeToMarkers[[cellType]] <- markers
}


markerset1 <- cellTypeToMarkers$Monocyte
markerset2 <- cellTypeToMarkers$Macrophage
markerset3 <- cellTypeToMarkers$Neutrophil
markerset4 <- cellTypeToMarkers$DC
combined_markers <- unique(c(markerset1, markerset2, markerset3, markerset4))

gene.list <- list(
    "B cell" = cellTypeToMarkers$'B cell',
    "Plasma cell" = cellTypeToMarkers$'Plasma cell',
    "T cell" = cellTypeToMarkers$'T cell',
    "Myeloid" = combined_markers,
    "NK" = cellTypeToMarkers$NK,
    "Mast" = cellTypeToMarkers$Mast
  )

seurat_immune_join <- JoinLayers(seurat_immune)
seurat_immune_join <- AddModuleScore(
    object = seurat_immune_join,
    features = gene.list,
    name = 'module_score_',
    layer = 'scale_data'
)

module_scores_df <- data.frame(
    ModuleScore1 = seurat_immune_join@meta.data$module_score_1,
    ModuleScore2 = seurat_immune_join@meta.data$module_score_2,
    ModuleScore3 = seurat_immune_join@meta.data$module_score_3,
    ModuleScore4 = seurat_immune_join@meta.data$module_score_4,
    ModuleScore5 = seurat_immune_join@meta.data$module_score_5,
    ModuleScore6 = seurat_immune_join@meta.data$module_score_6
 )

rownames(module_scores_df) <- colnames(seurat_immune_join)

clusters <- Idents(object = seurat_immune_join)
cell_ids_by_cluster <- split(x = colnames(seurat_immune_join), f = clusters)

average_scores_by_cluster <- list()
for (cluster in names(cell_ids_by_cluster)) {
  cell_ids <- cell_ids_by_cluster[[cluster]]
  cell_ids <- cell_ids[cell_ids %in% rownames(module_scores_df)]
  cluster_scores <- module_scores_df[cell_ids, , drop = FALSE]
  average_scores <- colMeans(cluster_scores, na.rm = TRUE)
  average_scores_by_cluster[[cluster]] <- average_scores
}

average_scores_df <- do.call(rbind, average_scores_by_cluster)
rownames(average_scores_df) <- names(average_scores_by_cluster)
colnames(average_scores_df) <- c("B cell", "Plasma cell", "T cell", "Myeloid", "NK", "Mast")
    
average_scores_df <- as.data.frame(scale(average_scores_df))

seurat_immune@meta.data$immune_subtype = ""
average_scores_df$type <- apply(average_scores_df, 1, function(row) {
    colnames(average_scores_df)[which.max(row)]
})

for(j in rownames(average_scores_df)) {
    type <- average_scores_df$type[rownames(average_scores_df) == j]
    seurat_immune@meta.data$immune_subtype[seurat_immune@meta.data$seurat_clusters == j] <- type}


pdf('immune_subtype_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_immune, reduction = "umap.mnn", label = TRUE, repel = TRUE,group.by='immune_subtype')
dev.off()

save(seurat_immune, file = 'seurat_immune.RData')

############################################################
#        extract the names of cells in each type
############################################################

B_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "B cell",])
Mast_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "Mast",])
Myeloid_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "Myeloid",])
NK_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "NK",])
Plasma_cell_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "Plasma cell",])
T_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "T cell",])


write(B_cells, file = "B cell.txt")
write(Mast_cells, file = "Mast.txt")
write(Myeloid_cells, file = "Myeloid.txt")
write(NK_cells, file = "NK.txt")
write(Plasma_cell_cells, file = "Plasma cell.txt")
write(T_cells, file = "T cells.txt")


ref <- readRDS(ref_immune.rds')
write.table(rownames(ref@meta.data[ref@meta.data$cell_type=='B cell',]),'B cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$cell_type=='mast cell',]),'Mast_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$cell_type=='myeloid cell',]),'Myeloid_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$cell_type=='natural killer cell',]),'NK_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$cell_type=='plasma cell',]),'Plasma cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$cell_type=='T cell',]),'T cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)


