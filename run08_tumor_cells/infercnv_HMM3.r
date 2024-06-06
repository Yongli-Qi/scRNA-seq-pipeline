library(infercnv, lib.loc = "/home/yq238/.conda/envs/scRNA/lib/R/library")
library(Seurat, lib.loc = "/home/yq238/.conda/envs/scRNA/lib/R/library")

options(scipen = 100)

setwd('/home/yq238/project/scRNA_seq/HTAN/08tumor_cells')
load("count_matrix_epithelial.RData")

col_names = colnames(count_matrix)

normal_cells <- col_names[grep("^HTA8_300", col_names)]
sclc_cells <- col_names[grep("^HTA8_2", col_names)]
nsclc_cells <- col_names[grep("^HTA8_1", col_names)]

anno_df <- data.frame(
  Cell_Name = col_names,
  Cell_Type = ifelse(col_names %in% normal_cells, "normal",
                     ifelse(col_names %in% sclc_cells, "sclc", 
                            ifelse(col_names %in% nsclc_cells, "nsclc", NA)))
)

rownames(anno_df) <- anno_df$Cell_Name
anno_df <- anno_df[-1]


gene_order <- read.csv('hg38_gencode_v27.csv',header = F)
rownames(gene_order)<-gene_order$V1
gene_order = gene_order[,2:4]

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = count_matrix,
                                    annotations_file = anno_df,
                                    gene_order_file = gene_order,
                                    ref_group_names = c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = "inferCNV_out_HMM3", 
                             cluster_by_groups=TRUE, 
                             plot_steps=FALSE,
                             denoise=TRUE,
                             HMM=TRUE,
			     HMM_type = 'i3',
                             no_prelim_plot=TRUE,
                             no_plot=TRUE,
                             num_threads = 5)  

save(infercnv_obj, file='infercnv_obj_HMM3.RData')
