library(Seurat)
library(dplyr)



#### calculate the number of cells after filtering empty drops
rdata_dir <- "/home/yq238/project/scRNA_seq/HTAN/01seurat_objs"
filtered_dir <- "/home/yq238/project/scRNA_seq/HTAN/02seurat_objs_filtered"
rdata_files <- list.files(rdata_dir, pattern = "\\.RData$", full.names = TRUE)

for (file in rdata_files) {
  load(file)
  sample <- sub("^.*seurat_object_(.*)\\.RData$", "\\1", basename(file))
  real_cells_file <- paste0(filtered_dir, "/", sample, "_real_cells.txt")
  real_cells <- read.table(real_cells_file, stringsAsFactors = FALSE, colClasses = "character")[, 1]
  seurat_obj_filtered <- subset(seurat_obj, cells = real_cells)
  cat(paste(basename(file), "- Cell count:", length(rownames(seurat_obj_filtered@meta.data)), "\n"))
}



#### calculate the number of cells after filtering based on "percent.mt"
data_dir <- "/home/yq238/project/scRNA_seq/HTAN2/02seurat_objs_filtered"

files <- list.files(path = data_dir, pattern = "\\.RData$", full.names = TRUE)
cell_counts <- numeric(length(files))

for (i in 1:length(files)) {
  load(files[i])
  cell_counts[i] <- length(rownames(seurat_obj_filtered@meta.data))
  cat(paste(basename(files[i]), "- Cell count:", cell_counts[i], "\n"))
}



























