setwd('/home/yq238/project/scRNA_seq/HTAN/06coarse_annotation')
load('/home/yq238/project/scRNA_seq/HTAN/05seurat_obj_integrated/seurat_integrated_fastmnn.RData')

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


total_variance <- seurat_integrated@reductions$pca@misc$total.variance
eigValues = (seurat_integrated[["pca"]]@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)

seurat_integrated_join <- JoinLayers(seurat_integrated)
save(seurat_integrated_join, file = 'seurat_integrated_fastmnn_join_2.6.RData')

pdf('gene_expression_umap_2.6.pdf', width = 7 , height = 6)
FeaturePlot(seurat_integrated_join, 
             reduction = "umap.mnn", 
             features = c("PTPRC","COL1A1","CLDN5","EPCAM"), 
             order = TRUE,
             min.cutoff = 'q10', 
             label = TRUE, raster=FALSE)
dev.off()


seurat_integrated@meta.data$coarse_annotation <- ifelse(seurat_integrated@meta.data$seurat_cluster %in% c(0,1,3,6,7,9,10,12,14,15,19,20,24,25,26,28,30,32,36,41,47,49,53,54,55,58,59,60), "Immune",
  ifelse(seurat_integrated@meta.data$seurat_cluster %in% c(17,27,46,61), "Mesenchymal", 
  ifelse(seurat_integrated@meta.data$seurat_cluster %in% c(2,4,5,8,11,13,16,18,21,22,23,29,31,33,34,35,37,38,39,40,42,43,44,45,48,50,51,52,56,57), "Epithelial", 
"NA")))

pdf('corse_annotation_umap_2.6.pdf', width = 7 , height = 6)
DimPlot(seurat_integrated, reduction = "umap.mnn", label = TRUE, repel = TRUE, group.by='coarse_annotation')
dev.off()

table(seurat_integrated@meta.data$coarse_annotation)
#Epithelial      Immune Mesenchymal 
#      78935       89191        7853 


#####################################################
      extract the names of cells in each type
#####################################################
cell_idents <- Idents(seurat_integrated_join)


immune_cells <- rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$coarse_annotation == "Immune",])
mesenchymal_cells <- rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$coarse_annotation == "Mesenchymal",])
epithelial_cells <- rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$coarse_annotation == "Epithelial",])

write(immune_cells, file = "immune_cells.txt")
write(mesenchymal_cells, file = "mesenchymal_cells.txt")
write(epithelial_cells, file = "epithelial_cells.txt")


save(seurat_integrated, file = 'seurat_integrated_coarse_annotation.RData')




