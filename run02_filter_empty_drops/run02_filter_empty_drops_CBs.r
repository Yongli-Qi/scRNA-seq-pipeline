library(Seurat)
library(scCB2)

setwd('/home/yq238/project/scRNA_seq/HTAN2')

##########################################
#         identify real cells
##########################################

input_dir <- "/home/yq238/project/scRNA_seq/HTAN/00data"
output_dir <- "/home/yq238/project/scRNA_seq/HTAN/02seurat_objs_filtered"
sub_dirs <- list.dirs(path = input_dir, full.names = TRUE, recursive = FALSE)
sub_dirs <- sub_dirs[-10]
sub_dirs <- sub_dirs[-10]
sub_dirs <- sub_dirs[-37]

for (dir in sub_dirs) {
  counts <- Read10X(data.dir = dir)
  counts <- as(counts, 'dgCMatrix')
  CBOut <- CB2FindCell(counts, FDR_threshold = 0.01,
                       lower = 100, Ncores = 2, verbose = TRUE)
  RealCell <- GetCellMat(CBOut)
  real_cells <- colnames(RealCell)
  sub_dir_name <- basename(dir)
  output_file <- paste0(output_dir, "/", sub_dir_name, "_real_cells.txt")
  write.table(data.frame(real_cells), file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#####################################################
#    filter seurat objects (CB2 and percent.mt 20)
#####################################################
rdata_dir <- "/home/yq238/project/scRNA_seq/HTAN/01seurat_objs"
filtered_dir <- "/home/yq238/project/scRNA_seq/HTAN/02seurat_objs_filtered"
rdata_files <- list.files(rdata_dir, pattern = "\\.RData$", full.names = TRUE)

for (file in rdata_files) {
  load(file)
  cat("Loaded:", file, "\n")
  sample <- sub("^.*seurat_object_(.*)\\.RData$", "\\1", basename(file))
  real_cells_file <- paste0(filtered_dir, "/", sample, "_real_cells.txt")
  real_cells <- read.table(real_cells_file, stringsAsFactors = FALSE, colClasses = "character")[, 1]
  seurat_obj_filtered <- subset(seurat_obj, cells = real_cells)
  seurat_obj_filtered[['percent.mt']] <- PercentageFeatureSet(seurat_obj_filtered,pattern = "^MT-")
  seurat_obj_filtered <- subset(x = seurat_obj_filtered, 
                         subset = (percent.mt < 20))
  save_file_name <- paste0(filtered_dir, "/seurat_object_", sample, "_filtered.RData")
  save(seurat_obj_filtered, file = save_file_name)
}







