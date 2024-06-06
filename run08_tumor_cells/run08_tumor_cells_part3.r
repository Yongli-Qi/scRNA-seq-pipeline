library(Seurat)

setwd('/home/yq238/project/scRNA_seq/HTAN/08tumor_cells')
load('infercnv_obj_HMM3.RData')

extract_cells_to_dataframe <- function(data) {
     results_df <- data.frame(Subcluster = character(), Cell_Name = character(), stringsAsFactors = FALSE)    

     for (subcluster in names(data)) {
         cell_names <- names(data[[subcluster]])
         temp_df <- data.frame(Subcluster = rep(subcluster, length(cell_names)), Cell_Name = cell_names, stringsAsFactors = FALSE)
         results_df <- rbind(results_df, temp_df)
     }
     
     return(results_df)
 }
 
nsclc_data <- infercnv_obj@tumor_subclusters$subclusters$nsclc
nsclc_dataframe <- extract_cells_to_dataframe(nsclc_data)
write.table(nsclc_dataframe, file='nsclc_cluster2cells.txt',quote = F,sep ='\t',col.names = F,row.names = F)

sclc_data <- infercnv_obj@tumor_subclusters$subclusters$sclc
sclc_dataframe <- extract_cells_to_dataframe(sclc_data)
write.table(sclc_dataframe, file='sclc_cluster2cells.txt',quote = F,sep ='\t',col.names = F,row.names = F)

normal_data <- infercnv_obj@tumor_subclusters$subclusters$normal
normal_dataframe <- extract_cells_to_dataframe(normal_data)
write.table(normal_dataframe, file='normal_cluster2cells.txt',quote = F,sep ='\t',col.names = F,row.names = F)
