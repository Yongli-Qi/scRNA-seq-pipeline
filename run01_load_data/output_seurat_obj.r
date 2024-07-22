library(argparse)
library(Seurat)


parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character")
args <- parser$parse_args()

create_seurat_object_from_dir <- function(dir) {
  model_name <- basename(dir)
  file1 <- file.path(dir, 'barcodes.tsv')
  file2 <- file.path(dir, 'genes.tsv')
  file3 <- file.path(dir, 'matrix.mtx')
  files_exist <- file.exists(file1, file2, file3)
  if (all(files_exist)) {
    counts <- Read10X(data.dir = dir)
    counts <- as.matrix(counts)
    seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 10, min.features = 100)
    save_file_name <- paste0("seurat_object_", model_name, ".RData")
    save(seurat_obj, file = save_file_name)
  } else {
    message("One or more files do not exist in the specified directory.")
  }
}

create_seurat_object_from_dir(args$input)