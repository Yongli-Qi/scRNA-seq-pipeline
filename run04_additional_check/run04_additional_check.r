library(Seurat)

################################################
#    creat a list for storing seurat objs
################################################
setwd("/home/yq238/project/scRNA_seq/HTAN/03seurat_objs_removed")
rdata_files <- list.files(pattern = "\\.RData$")
seurat_objects <- list()
for (file in rdata_files) {
  model_name <- sub("^.*seurat_object_(.*)\\_filtered.RData$", "\\1", basename(file))
  load(file)
  cat("Loaded:", file, "\n")
  seurat_objects[[model_name]] <- seurat_obj_removed
  }


setwd("/home/yq238/project/scRNA_seq/HTAN/04additional_check")
samples <- c('HTA8_1001_1','HTA8_1002_1','HTA8_1003_1','HTA8_1004_1','HTA8_1005_1','HTA8_1005_2','HTA8_1006_1','HTA8_1007_1','HTA8_1008_1','HTA8_1009_1','HTA8_1009_2','HTA8_1009_3','HTA8_1010_1','HTA8_1011_1','HTA8_1012_1','HTA8_1013_1','HTA8_1014_1','HTA8_1015_1','HTA8_1016_1','HTA8_1017_1','HTA8_1018_1','HTA8_1019_1','HTA8_1020_1','HTA8_1021_1','HTA8_1022_1','HTA8_2001_1','HTA8_2002_1','HTA8_2003_1','HTA8_2004_1','HTA8_2005_1','HTA8_2005_2','HTA8_2005_3','HTA8_2006_1','HTA8_2007_1','HTA8_2007_2','HTA8_2008_1','HTA8_2009_1','HTA8_2010_1','HTA8_2010_2','HTA8_2011_1','HTA8_2012_1','HTA8_2013_1','HTA8_2014_1','HTA8_2015_1','HTA8_2016_1','HTA8_2017_1','HTA8_2018_1','HTA8_2019_1','HTA8_3001_1','HTA8_3002_1','HTA8_3003_1','HTA8_3004_1')
seurat_combined = merge(seurat_objects[[1]],y = seurat_objects[-1],add.cell.ids = samples)

orig <- NULL
for(i in 1:length(seurat_objects)){
  orig = append(orig,rep(samples[i],ncol(seurat_objects[[i]])))
}
seurat_combined[["orig.ident"]] = orig

seurat_combined[['percent.mt']] <- PercentageFeatureSet(seurat_combined,pattern = "^MT-")

save(seurat_combined, file = 'seurat_combined.RData')


################################################
#             check the distribution
################################################
new_folder_name <- "/home/yq238/project/scRNA_seq/HTAN/04additional_check"
seurat_objects <- setNames(
  lapply(names(seurat_objects), function(name) {
    x <- seurat_objects[[name]]
    plot <- VlnPlot(x, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
    plot_filename <- paste0(new_folder_name, "/", name, "_plot.png")
    ggsave(plot_filename, plot = plot, width = 10, height = 6)
    return(x)
  }), 
  names(seurat_objects)
)


################################################
#      normalization, scaling, and PCA
################################################
seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_combined <- FindVariableFeatures(seurat_combined)
seurat_combined <- ScaleData(seurat_combined)
seurat_combined <- RunPCA(seurat_combined)

save(seurat_combined, file = 'seurat_combined_PCA.RData')

# plot the most variable genes
Top10 <- head(VariableFeatures(seurat_combined), 10)
pdf('variable_genes.pdf', width = 7 , height = 6)
plot1 <- VariableFeaturePlot(seurat_combined)
plot2 <- LabelPoints(plot=plot1, points = Top10, repel = TRUE,xnudge = 0,ynudge = 0)
plot2
dev.off()

################################################
#    check cell cycle phase and percent.mt
################################################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_combined <- JoinLayers(seurat_combined)
seurat_combined <- CellCycleScoring(seurat_combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Plot the PCA colored by cell cycle phase
pdf('pca_cell_cycle_phase.pdf', width = 7 , height = 6)
DimPlot(seurat_combined,reduction = 'pca',group.by = 'Phase',split.by = 'Phase')
dev.off()
# Create a new column, dividing percent.mt into two groups
seurat_combined$mt_group <- ifelse(seurat_combined$percent.mt > median(seurat_combined$percent.mt), "High_MT", "Low_MT")
# Plot a PCA plot, coloring by mt_group
pdf('pca_mt_group.pdf', width = 7 , height = 6)
DimPlot(seurat_combined,reduction = 'pca',group.by = 'mt_group',split.by = 'mt_group')
dev.off()



