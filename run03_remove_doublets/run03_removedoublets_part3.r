############################################
#   R: subset the non-doublet cells
##############################################
library(Seurat)

setwd("/home/yq238/project/scRNA_seq/HTAN/03seurat_objs_removed")

rdata_dir <- "/home/yq238/project/scRNA_seq/HTAN/02seurat_objs_filtered"
filtered_dir <- "/home/yq238/project/scRNA_seq/HTAN/03seurat_objs_removed"
rdata_files <- list.files(rdata_dir, pattern = "\\.RData$", full.names = TRUE)

for (file in rdata_files) {
  load(file)
  cat("Loaded:", file, "\n")
  sample <- sub("^.*seurat_object_(.*)\\_filtered.RData$", "\\1", basename(file))
  non_doublets_cells_file <- paste0(filtered_dir, "/", sample, "_non_doublets_cells.txt")
  non_doublets_cells <- read.table(non_doublets_cells_file, stringsAsFactors = FALSE, colClasses = "character")[, 1]
  seurat_obj_removed <- subset(seurat_obj_filtered, cells = non_doublets_cells)
  save_file_name <- paste0(filtered_dir, "/seurat_object_", sample, "_removed.RData")
  save(seurat_obj_removed, file = save_file_name)
}







