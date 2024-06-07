data_dir <- "/home/yq238/project/scRNA_seq/HTAN/03seurat_objs_removed"

files <- list.files(path = data_dir, pattern = "\\.RData$", full.names = TRUE)
cell_counts <- numeric(length(files))

for (i in 1:length(files)) {
  load(files[i])
  cell_counts[i] <- length(rownames(seurat_obj_removed@meta.data))
  cat(paste(basename(files[i]), "- Cell count:", cell_counts[i], "\n"))
}



























