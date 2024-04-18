######################################################
# extract count matrix in T cells
######################################################
library(Seurat)
library(SeuratWrappers)

setwd('/home/yq238/project/scRNA_seq/HTAN2/08T_cell_subtype_revised')
load('/home/yq238/project/scRNA_seq/HTAN2/07immune_subtype_revised/seurat_immune.RData')

seurat_T <- subset(x = seurat_immune, subset = (immune_subtype == 'T cell'))
seurat_T_join <- JoinLayers(seurat_T)
count_matrix <- GetAssayData(seurat_T_join, layer = "counts")
#count_matrix <- as.matrix(count_matrix)
#write.table(count_matrix, "count_matrix_unnormalized_T.txt", sep ='\t', quote =FALSE)
save(count_matrix, file='count_matrix_unnormalized_T_sparse.RData')

######################################################
# perform NMF method using RcppML
######################################################
library(RcppML)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggrepel)
library(uwot)
library(dplyr)

load('count_matrix_unnormalized_T_sparse.RData')

# decide how many factors we are going to use
crossValidate_result <- crossValidate(count_matrix, k = c(1,5,10,15,20,25,30,35,40,45,50), reps = 1, verbose = TRUE)
crossValidate_result2 <- crossValidate(count_matrix, k = c(2,3,4,21,22,23,24,26,27,28,29), reps = 1, verbose = TRUE)
load('crossValidate_result_part1.RData')
load('crossValidate_result_part2.RData')

crossValidate_result <- cbind(crossValidate_result, crossValidate_result2)
plot(crossValidate_result)

model <- RcppML::nmf(count_matrix, k = 30, maxit = 500)

gene_loadings <- model@w

calculate_sums <- function(gene_set, gene_loadings) {
  subset_loadings <- gene_loadings[rownames(gene_loadings) %in% gene_set, ]
  colSums(subset_loadings)
}

Exhausted_markers <- c('CTLA4','LAG3','PTMS','PDCD1','TIGIT')
Effector_markers <- c('CD8A','CCL4','CST7','GZMA','GZMH','CTSW','PRF1','GZMB','FGFBP2','FCGR3A','TRGC2','GNLY','TYROBP')
Memory_markers <- c('CD44','GZMK','IL7R','CD69','CD27')
Tgd_markers <-c('FCER1G','CMC1','KLRD1','KLRG1','KLRF1','TRDC')
Treg_markers <- c('FOXP3','TNFRSF4','CARD16','IL2RA','TBC1D4')
Tconv_markers <- c('CD4','LTB','SPOCK2','SOCS3','PBXIP1')


Tgd_set <- calculate_sums(Tgd_markers, gene_loadings)
Effector_set <- calculate_sums(Effector_markers, gene_loadings)
Memory_set <- calculate_sums(Memory_markers, gene_loadings)
Exhausted_set <- calculate_sums(Exhausted_markers, gene_loadings)
Treg_set <- calculate_sums(Treg_markers, gene_loadings)
Tconv_set <- calculate_sums(Tconv_markers, gene_loadings)


total_loadings <- data.frame(Tgd = Tgd_set, Effector = Effector_set, Memory = Memory_set, Exhausted = Exhausted_set, Treg = Treg_set, Tconv = Tconv_set)

cell_loadings <- model@h
selected_nmf <- cell_loadings[c("nmf22", "nmf21", "nmf9", "nmf7", "nmf8", "nmf2","nmf12"),]

max_nmf_df <- data.frame(cell = colnames(selected_nmf), nmf = NA, stringsAsFactors = FALSE)

for (i in 1:ncol(selected_nmf)) {
  max_nmf <- which.max(selected_nmf[, i])
  max_nmf_df$nmf[i] <- rownames(selected_nmf)[max_nmf]
}

max_nmf_df$T_subtype <- ifelse(max_nmf_df$nmf == "nmf7", "CD8+ T Exhausted",
                               ifelse(max_nmf_df$nmf == "nmf9", "CD8+ T Memory-like",
                                      ifelse(max_nmf_df$nmf == "nmf21", "CD8+ T Effector-like",
                                             ifelse(max_nmf_df$nmf == "nmf22", "Tgd",
                                                    ifelse(max_nmf_df$nmf == "nmf8", "CD4+ Treg",
                                                           ifelse(max_nmf_df$nmf == "nmf2", "CD4+ Tconv", 
                                                                  ifelse(max_nmf_df$nmf == "nmf12", "CD4+ Tconv",
                                                                  NA)))))))
write.table(max_nmf_df,file='T_subtype.txt',sep = '\t',quote = F,row.names = F)


######################################################
# extract cell names from reference
######################################################


ref <- readRDS('T_ref.rds')
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='Tgd',]),'Tgd_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='CD8+ T Effector-like',]),'Effector_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)

write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='CD8+ T Memory-like',]),'Memory_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)

write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='CD8+ T Exhausted',]),'Exhausted_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='CD4+ Treg',]),'Treg_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)
write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='CD4+ Tconv',]),'Tconv_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)


######################################################
# add annotation in seurat obj and visulize
######################################################

load("/home/yq238/project/scRNA_seq/HTAN2/07immune_subtype_revised/seurat_immune.RData")
seurat_T <- subset(x = seurat_immune, subset = (immune_subtype == 'T cell'))

T_subtype <- read.table("T_subtype.txt",sep = '\t',header = T)
Tgd <- T_subtype$cell[T_subtype$T_subtype=='Tgd']
Effector <- T_subtype$cell[T_subtype$T_subtype=='CD8+ T Effector-like']
Memory <- T_subtype$cell[T_subtype$T_subtype=='CD8+ T Memory-like']
Exhausted <- T_subtype$cell[T_subtype$T_subtype=='CD8+ T Exhausted']
Treg <- T_subtype$cell[T_subtype$T_subtype=='CD4+ Treg']
Tconv <- T_subtype$cell[T_subtype$T_subtype=='CD4+ Tconv']


seurat_T@meta.data$T_subtype <- ifelse(rownames(seurat_T@meta.data) %in% Tgd, "Tgd",
ifelse(rownames(seurat_T@meta.data) %in% Effector, "CD8+ T Effector-like",
ifelse(rownames(seurat_T@meta.data) %in% Memory, "CD8+ T Memory-like",
ifelse(rownames(seurat_T@meta.data) %in% Exhausted, "CD8+ T Exhausted",
ifelse(rownames(seurat_T@meta.data) %in% Treg, "CD4+ Treg",
ifelse(rownames(seurat_T@meta.data) %in% Tconv, "CD4+ Tconv",NA
))))))

save(seurat_T, file="seurat_T_subtype_annotation.RData")


seurat_T <- FindVariableFeatures(seurat_T)
seurat_T <- ScaleData(seurat_T)
seurat_T <- RunPCA(seurat_T)
seurat_T <- IntegrateLayers(object = seurat_T, method = FastMNNIntegration, new.reduction = "integrated.mnn", verbose = FALSE)

seurat_T <- FindNeighbors(seurat_T, reduction = "integrated.mnn", dims = 1:50)
seurat_T <- FindClusters(seurat_T, resolution = 0.8, cluster.name = "mnn_clusters")
seurat_T <- RunUMAP(seurat_T, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")


pdf('T_cells_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_T, reduction = "umap.mnn", label = TRUE, repel = TRUE, group.by='T_subtype')
dev.off()

markers <- read.table('markers_T.txt')
markers <- markers$V1

seurat_T_join <- JoinLayers(seurat_T)

desired_order <- c("Tgd","CD8+ T Effector-like", "CD8+ T Memory-like","CD8+ T Exhausted","CD4+ Treg", "CD4+ Tconv")
seurat_T_join$T_subtype <- factor(seurat_T_join$T_subtype, levels = desired_order)

pdf('T_cells_expression_dotplot_scaled.pdf', width = 14 , height = 6)
DotPlot(seurat_T_join, features = markers, assay = 'RNA', group.by = 'T_subtype', cols = c('white','red')) + RotatedAxis()
dev.off()

save(seurat_T_join, file='seurat_T_join_forplot.RData')



