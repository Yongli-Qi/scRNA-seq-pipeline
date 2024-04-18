library(Seurat)
library(ggplot2)

setwd('/home/yq238/project/scRNA_seq/HTAN2/15heatmap_SCLC')

load('/home/yq238/project/scRNA_seq/HTAN2/15DEG_SCLC/seurat_SCLC_join.RData')


tf <- read.table('TFs.txt')
tf <- tf$V1

seurat_SCLC_join <- ScaleData(seurat_SCLC_join, features = rownames(seurat_SCLC_join))

pdf('heatmap_TFs_scale.pdf',width= 7, height = 12)
tf_heatmap_scale <- DoHeatmap(seurat_SCLC_join, features = tf, group.by = "ident")
tf_heatmap_scale
dev.off()

pdf('heatmap_TFs_scale_nogene.pdf',width= 7, height = 12)
tf_heatmap_scale + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

write.table(tf_heatmap_scale[[1]][["data"]], file = "tf_heatmap_scaled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)

#############################################################################

er <- read.table('EpiGenes_main.txt')
er <- er$V1

pdf('heatmap_ERs_scale_nogene.pdf',width= 7, height = 12)
er_heatmap_scale <- DoHeatmap(seurat_SCLC_join, features = er, group.by = "ident")
er_heatmap_scale + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

write.table(er_heatmap_scale[[1]][["data"]], file = "er_heatmap_scaled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)

############################################################################
genes <- c('HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'TAP2', 'SOX2', 'MYC', 'MYCL', 'MYCN', 'MAX')

pdf('heatmap_seleted_geness_scale.pdf',width= 7, height = 12)
genes_heatmap_scale <- DoHeatmap(seurat_SCLC_join, features = genes, group.by = "ident")
dev.off()

write.table(genes_heatmap_scale[[1]][["data"]], file = "seleted_genes_heatmap_scaled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)



