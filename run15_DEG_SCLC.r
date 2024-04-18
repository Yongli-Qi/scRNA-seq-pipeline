setwd('/home/yq238/project/scRNA_seq/HTAN2/15DEG_SCLC')

load('/home/yq238/project/scRNA_seq/HTAN2/13cancer_cells_violin/seurat_SCLC_join_rePCA.RData')
seurat_SCLC_join <- FindNeighbors(seurat_SCLC_join, reduction = "pca", dims = 1:50)
seurat_SCLC_join <- FindClusters(seurat_SCLC_join, resolution = 0.4, cluster.name = "clusters")
seurat_SCLC_join <- RunUMAP(seurat_SCLC_join, reduction = "pca", dims = 1:50, reduction.name = "umap")
seurat_SCLC_join@meta.data$SCLC_subtype <- ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(0,5,6,7,8,9,10,11,13,14,16,20,21,22), "SCLC-A",
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(2,3,4,15), "SCLC-N", 
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(12), "SCLC-P", 
  ifelse(seurat_SCLC_join@meta.data$seurat_cluster %in% c(18,19), "SCLC-Y",
"NA"))))

Idents(seurat_SCLC_join) <- seurat_SCLC_join@meta.data$SCLC_subtype

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-A', ident.2 = 'SCLC-N',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_A_vs_N.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-A', ident.2 = 'SCLC-P',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_A_vs_P.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-A', ident.2 = 'SCLC-Y',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_A_vs_Y.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-N', ident.2 = 'SCLC-P',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_N_vs_P.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-N', ident.2 = 'SCLC-Y',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_N_vs_Y.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-P', ident.2 = 'SCLC-Y',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_P_vs_Y.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'NA', assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_NA_vs_rest.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-A', assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_A_vs_rest.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-N', assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_N_vs_rest.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-P', assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_P_vs_rest.csv',quote = F)

deg <- FindMarkers(seurat_SCLC_join, ident.1 = 'SCLC-Y', assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_Y_vs_rest.csv',quote = F)


####################################################
load("../10_epithelial_tumors_raw_counts/seurat_epithelial_tumorstatus_annotation.RData")
seurat_epithelial_join <- JoinLayers(seurat_epithelial)

seurat_tumor <- subset(x = seurat_epithelial_join, subset = (tumorstatus == 'cancer'))

seurat_tumor@meta.data$cancer_type <- ifelse(seurat_tumor@meta.data$ne_subtype == "non-neuroendocrine", "NSCLC", "SCLC")

Idents(seurat_tumor) <- seurat_tumor@meta.data$cancer_type

deg <- FindMarkers(seurat_tumor, ident.1 = 'SCLC', ident.2 = 'NSCLC',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_vs_NSCLC.csv',quote = F)

####################################################
load("/gpfs/gibbs/project/augert/yq238/scRNA_seq/HTAN2/06coarse_annotation/seurat_integrated_coarse_annotation.RData")
load('/home/yq238/project/scRNA_seq/HTAN2/07immune_subtype_revised/seurat_immune.RData')
load('/home/yq238/project/scRNA_seq/HTAN2/08T_cell_subtype_revised/seurat_T_subtype_annotation.RData')
load('/home/yq238/project/scRNA_seq/HTAN2/10_epithelial_tumors_raw_counts/seurat_epithelial_tumorstatus_annotation.RData')
load('/home/yq238/project/scRNA_seq/HTAN2/11_normal_non_ne_subtype_revised_allnonnecells/seurat_epithelial_non_ne_normal_subtype.RData')
load('/home/yq238/project/scRNA_seq/HTAN2/12extract_NSCLC_SCLC/seurat_NSCLC_join.RData')
load('/home/yq238/project/scRNA_seq/HTAN2/13cancer_cells_violin/seurat_SCLC_join_rePCA.RData')

B_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "B cell",])
Mast_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "Mast",])
Myeloid_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "Myeloid",])
NK_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "NK",])
Plasma_cells <- rownames(seurat_immune@meta.data[seurat_immune@meta.data$immune_subtype == "Plasma cell",])

Tgd_cells <- rownames(seurat_T@meta.data[seurat_T@meta.data$T_subtype == "Tgd",])
Effector_cells <- rownames(seurat_T@meta.data[seurat_T@meta.data$T_subtype == "CD8+ T Effector-like",])
Memory_cells <- rownames(seurat_T@meta.data[seurat_T@meta.data$T_subtype == "CD8+ T Memory-like",])
Exhausted_cells <- rownames(seurat_T@meta.data[seurat_T@meta.data$T_subtype == "CD8+ T Exhausted",])
Treg_cells <- rownames(seurat_T@meta.data[seurat_T@meta.data$T_subtype == "CD4+ Treg",])
Tconv_cells <- rownames(seurat_T@meta.data[seurat_T@meta.data$T_subtype == "CD4+ Tconv",])

Mesenchymal_cells <- rownames(seurat_integrated@meta.data[seurat_integrated@meta.data$coarse_annotation == "Mesenchymal",])
NSCLC_cells <- rownames(seurat_NSCLC_join@meta.data)

SCLC_cells <- rownames(seurat_SCLC_join@meta.data)
AE1_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "AE1",])
AE2_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "AE2",])
Basal_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Basal",])
Ciliated_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Ciliated",])
Club_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Club",])
Hepatocyte_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Hepatocyte",])
Ionocyte_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Ionocyte",])
Mucinous_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Mucinous",])
Tuft_cells <- rownames(seurat_epithelial_non_ne_normal@meta.data[seurat_epithelial_non_ne_normal@meta.data$non_ne_subtype == "Tuft",])

seurat_epithelial_normal <- subset(x = seurat_epithelial, subset = (tumorstatus == 'normal'))
seurat_epithelial_normal_ne <- subset(x = seurat_epithelial_normal, subset = (ne_subtype == 'neuroendocrine'))
neuroendocrine_cells <- rownames(seurat_epithelial_normal_ne@meta.data)

seurat_integrated@meta.data$final_annotated_type <- NA

seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% B_cells] <- 'B cell'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Mast_cells] <- 'Mast'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Myeloid_cells] <- 'Myeloid'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% NK_cells] <- 'NK'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Plasma_cells] <- 'Plasma cell'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Tgd_cells] <- 'Tgd'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Effector_cells] <- 'CD8+ T Effector-like'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Memory_cells] <- 'CD8+ T Memory-like'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Exhausted_cells] <- 'CD8+ T Exhausted'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Treg_cells] <- 'CD4+ Treg'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Tconv_cells] <- 'CD4+ Tconv'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Mesenchymal_cells] <- 'Mesenchymal'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% NSCLC_cells] <- 'NSCLC'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% SCLC_cells] <- 'SCLC'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% AE1_cells] <- 'AE1'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% AE2_cells] <- 'AE2'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Basal_cells] <- 'Basal'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Ciliated_cells] <- 'Ciliated'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Club_cells] <- 'Club'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Hepatocyte_cells] <- 'Hepatocyte'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Ionocyte_cells] <- 'Ionocyte'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Mucinous_cells] <- 'Mucinous'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% Tuft_cells] <- 'Tuft'
seurat_integrated@meta.data$final_annotated_type[rownames(seurat_integrated@meta.data) %in% neuroendocrine_cells] <- 'Neuroendocrine'


seurat_integrated <- JoinLayers(seurat_integrated)
Idents(seurat_integrated) <- seurat_integrated@meta.data$final_annotated_type

deg <- FindMarkers(seurat_integrated, ident.1 = 'SCLC', ident.2 = 'Neuroendocrine',  assay = "RNA", layer = 'data')
write.csv(deg, file = 'SCLC_vs_Neuroendocrine.csv',quote = F)



















