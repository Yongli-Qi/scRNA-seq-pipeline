library(Seurat)
.libPaths(c("/vast/palmer/apps/avx2/software/R/4.3.0-foss-2022b/lib64/R/library", .libPaths()))

setwd('/home/yq238/project/scRNA_seq/HTAN2/11_normal_non_ne_subtype_revised_allnonnecells')
load("../10_epithelial_tumors_raw_counts/seurat_epithelial_tumorstatus_annotation.RData")

seurat_epithelial_non_ne <- subset(x = seurat_epithelial, subset = (ne_subtype == 'non-neuroendocrine'))

seurat_epithelial_non_ne <- FindVariableFeatures(seurat_epithelial_non_ne)
seurat_epithelial_non_ne <- ScaleData(seurat_epithelial_non_ne)
seurat_epithelial_non_ne <- RunPCA(seurat_epithelial_non_ne)
seurat_epithelial_non_ne <- IntegrateLayers(object = seurat_epithelial_non_ne, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)

save(seurat_epithelial_non_ne, file="seurat_epithelial_non_ne_integrated.RData")

seurat_epithelial_non_ne <- FindNeighbors(seurat_epithelial_non_ne, reduction = "integrated.mnn", dims = 1:50)
seurat_epithelial_non_ne <- FindClusters(seurat_epithelial_non_ne, resolution = 3, cluster.name = "mnn_clusters")
seurat_epithelial_non_ne <- RunUMAP(seurat_epithelial_non_ne, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")


pdf('dimplot_clusters_epithelial_non_ne_mnn_3.pdf', width = 7 , height = 6)
DimPlot(
  seurat_epithelial_non_ne,
  reduction = "umap.mnn",
  group.by = c("mnn_clusters"),
  combine = FALSE, label.size = 2
)
dev.off()


total_variance <- seurat_epithelial_non_ne@reductions$pca@misc$total.variance
eigValues = (seurat_epithelial_non_ne[["pca"]]@stdev)^2
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

gene.list <- list(
    "AE1" = cellTypeToMarkers$AE1,
    "AE2" = cellTypeToMarkers$AE2,
    "Basal" = cellTypeToMarkers$Basal,
    "Ciliated" = cellTypeToMarkers$Ciliated,
    "Club" = cellTypeToMarkers$Club,
    "Hepatocyte" = cellTypeToMarkers$Hepatocyte,
    "Ionocyte" = cellTypeToMarkers$Ionocyte,
    "Mucinous" = cellTypeToMarkers$Mucinous,
    "Tuft" = cellTypeToMarkers$Tuft
)

seurat_epithelial_non_ne_join <- JoinLayers(seurat_epithelial_non_ne)
seurat_epithelial_non_ne_join <- AddModuleScore(
    object = seurat_epithelial_non_ne_join,
    features = gene.list,
    name = 'module_score_',
    layer = 'scale_data'
)


module_scores_df <- data.frame(
    ModuleScore1 = seurat_epithelial_non_ne_join@meta.data$module_score_1,
    ModuleScore2 = seurat_epithelial_non_ne_join@meta.data$module_score_2,
    ModuleScore3 = seurat_epithelial_non_ne_join@meta.data$module_score_3,
    ModuleScore4 = seurat_epithelial_non_ne_join@meta.data$module_score_4,
    ModuleScore5 = seurat_epithelial_non_ne_join@meta.data$module_score_5,
    ModuleScore6 = seurat_epithelial_non_ne_join@meta.data$module_score_6,
    ModuleScore7 = seurat_epithelial_non_ne_join@meta.data$module_score_7,
    ModuleScore8 = seurat_epithelial_non_ne_join@meta.data$module_score_8,
    ModuleScore9 = seurat_epithelial_non_ne_join@meta.data$module_score_9
 )

rownames(module_scores_df) <- colnames(seurat_epithelial_non_ne_join)

clusters <- Idents(object = seurat_epithelial_non_ne_join)
cell_ids_by_cluster <- split(x = colnames(seurat_epithelial_non_ne_join), f = clusters)

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
colnames(average_scores_df) <- c("AE1", "AE2", "Basal", "Ciliated", "Club", "Hepatocyte","Ionocyte","Mucinous","Tuft")

average_scores_df <- as.data.frame(scale(average_scores_df))

seurat_epithelial_non_ne@meta.data$non_ne_subtype = ""
average_scores_df$type <- apply(average_scores_df, 1, function(row) {
    colnames(average_scores_df)[which.max(row)]
})

for(j in rownames(average_scores_df)) {
    type <- average_scores_df$type[rownames(average_scores_df) == j]
    seurat_epithelial_non_ne@meta.data$non_ne_subtype[seurat_epithelial_non_ne@meta.data$seurat_clusters == j] <- type}


pdf('epithelial_non_ne_subtype_annotation_umap_2.pdf', width = 7 , height = 6)
DimPlot(seurat_epithelial_non_ne, reduction = "umap.mnn", label = TRUE, repel = TRUE,group.by='non_ne_subtype')
dev.off()


seurat_epithelial_non_ne_normal <- subset(x = seurat_epithelial_non_ne, subset = (tumorstatus == 'normal'))
pdf('epithelial_non_ne_normal_subtype_annotation_umap_2.pdf', width = 7 , height = 6)
DimPlot(seurat_epithelial_non_ne_normal, reduction = "umap.mnn", label = TRUE, repel = TRUE,group.by='non_ne_subtype')
dev.off()


table(seurat_epithelial_non_ne_normal@meta.data$non_ne_normal_subtype)

#3
#       AE1        AE2      Basal   Ciliated      Club    Hepatocyte   Ionocyte     Mucinous       Tuft 
#       549        734       7083        563        423       8161       1103         2444        375 



##############extract the cells

AE1_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "AE1",])
write(AE1_cells, file = "AE1 cell.txt")

AE2_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "AE2",])
write(AE2_cells, file = "AE2 cell.txt")

Basal_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Basal",])
write(Basal_cells, file = "Basal cell.txt")

Ciliated_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Ciliated",])
write(Ciliated_cells, file = "Ciliated cell.txt")

Club_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Club",])
write(Club_cells, file = "Club cell.txt")

Hepatocyte_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Hepatocyte",])
write(Hepatocyte_cells, file = "Hepatocyte cell.txt")

Ionocyte_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Ionocyte",])
write(Ionocyte_cells, file = "Ionocyte cell.txt")

Mucinous_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Mucinous",])
write(Mucinous_cells, file = "Mucinous cell.txt")

Tuft_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Tuft",])
write(Tuft_cells, file = "Tuft cell.txt")


##############extract the cells in ref

ref <- readRDS('/home/yq238/project/scRNA_seq/HTAN2/09_epithelial_subtype/ref_epi.rds')

write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='AE1',]),'AE1_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Ionocyte',]),'Ionocyte_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Hepatocyte',]),'Hepatocyte_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Mucinous',]),'Mucinous_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='AE2',]),'AE2_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Basal',]),'Basal_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Ciliated',]),'Ciliated_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Club',]),'Club_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Tuft',]),'Tuft_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)


save(seurat_epithelial_non_ne_normal, file="seurat_epithelial_non_ne_normal_subtype.RData")
















