############################################
#   extract the loaded cells
##############################################
library(Seurat)
setwd("/home/yq238/project/scRNA_seq/HTAN/02seurat_objs_filtered")
rdata_dir <- "/home/yq238/project/scRNA_seq/HTAN/02seurat_objs_filtered"
rdata_files <- list.files(rdata_dir, pattern = "\\.RData$", full.names = TRUE)

for (file in rdata_files) {
  load(file)
  cat("Loaded:", file, "\n")
  sample <- sub("^.*seurat_object_(.*)\\_filtered.RData$", "\\1", basename(file))
  save_file_name <- paste0(sample, "_cells_after_filtering.txt")
  write.table(rownames(seurat_obj_filtered@meta.data),save_file_name,sep = '\t',quote = F, col.names = F, row.names = F)
}







