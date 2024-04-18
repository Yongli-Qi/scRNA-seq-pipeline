library(Seurat)

setwd('/home/yq238/project/scRNA_seq/HTAN2/10_epithelial_tumors_raw_counts')
load("../09_epithelial_subtype_revised/seurat_epithelial_integrated.RData")

seurat_epithelial_join <- JoinLayers(seurat_epithelial)
count_matrix <- GetAssayData(seurat_epithelial_join, layer = "counts")
save(count_matrix, file="count_matrix_epithelial_sparse.RData")

#count_matrix <- as.matrix(count_matrix)
#save(count_matrix, file="count_matrix_epithelial.RData")


########################################
********* run infercnv HMM3 *************
########################################

############ R script to perform the analysis

library(infercnv, lib.loc = "/home/yq238/.conda/envs/scRNA/lib/R/library")
library(Seurat, lib.loc = "/home/yq238/.conda/envs/scRNA/lib/R/library")

options(scipen = 100)

setwd('/home/yq238/project/scRNA_seq/HTAN2/10_epithelial_tumors_raw_counts')
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


############ bash script to submit the job

#!/bin/bash
#SBATCH --job-name=infercnv_HMM3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yongli.qi@yale.edu
#SBATCH -p day
#SBATCH --cpus-per-task=5
#SBATCH --mem=500G
#SBATCH -t 24:00:00

module load R/4.3.0-foss-2022b

export PATH=/home/yq238/miniconda3/bin:$PATH
source activate scRNA

export LD_LIBRARY_PATH=/home/yq238/.conda/envs/scRNA/lib/JAGS/modules-4:$LD_LIBRARY_PATH

Rscript infercnv_HMM3.r


######################################################
# extract cells in each cluster identified by infercnv
######################################################

load('infercnv_obj_HMM3.RData')

nsclc_data <- infercnv_obj@tumor_subclusters$subclusters$nsclc
extract_cells_to_dataframe <- function(nsclc_data) {
     results_df <- data.frame(Subcluster = character(), Cell_Name = character(), stringsAsFactors = FALSE)    

     for (subcluster in names(nsclc_data)) {
         cell_names <- names(nsclc_data[[subcluster]])
         temp_df <- data.frame(Subcluster = rep(subcluster, length(cell_names)), Cell_Name = cell_names, stringsAsFactors = FALSE)
         results_df <- rbind(results_df, temp_df)
     }
     
     return(results_df)
 }
 
nsclc_dataframe <- extract_cells_to_dataframe(nsclc_data)
write.table(nsclc_dataframe, file='nsclc_cluster2cells.txt',quote = F,sep ='\t',col.names = F,row.names = F)

sclc_data <- infercnv_obj@tumor_subclusters$subclusters$sclc
sclc_dataframe <- extract_cells_to_dataframe(sclc_data)
write.table(sclc_dataframe, file='sclc_cluster2cells.txt',quote = F,sep ='\t',col.names = F,row.names = F)

normal_data <- infercnv_obj@tumor_subclusters$subclusters$normal
normal_dataframe <- extract_cells_to_dataframe(normal_data)
write.table(normal_dataframe, file='normal_cluster2cells.txt',quote = F,sep ='\t',col.names = F,row.names = F)


######################################################
# extract cnv percentage of cells using python
######################################################

########## Only genes identified in cBioPortal as having
########## mutations in SCLC and NSCLC are retained in our analysis

a0=open('shared_genes.txt','r')
a1=open('HMM_CNV_predictions.HMMi3.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat','r')
w=open('HMM3_CNV_filtered.txt','w')

genes_list=[]
for line in a0:
    li=line.split()
    genes_list.append(li[0])

for line in a1:
    li=line.split()
    if li[3] in genes_list:
        w.write(line)

w.close()

##########extract mutated cells

def extract_cluster2cells(file_in):
    a=open(file_in,'r')

    dict_cluster2cells={}
    for line in a:
        li=line.split()
        cluster=li[0]
        cell=li[1]
        if cluster not in dict_cluster2cells:
            list=[]
            list.append(li[1])
            dict_cluster2cells[cluster] = list
        else:
            dict_cluster2cells[cluster].append(li[1])
    return dict_cluster2cells

dict_normal = extract_cluster2cells('normal_cluster2cells.txt')
dict_sclc = extract_cluster2cells('sclc_cluster2cells.txt')
dict_nsclc = extract_cluster2cells('nsclc_cluster2cells.txt')


a1=open('HMM3_CNV_filtered.txt','r')
w1=open('cells_cnv_sclc.txt','w')
w2=open('cells_cnv_nsclc.txt','w')
w3=open('cells_cnv_normal.txt','w')

normal_cells=[]
sclc_cells=[]
nsclc_cells=[]

for line in a1:
    li=line.split()
    cluster = li[0].split('.')[1]
    if cluster.startswith('sclc'):
        for cell in dict_sclc[cluster]:
            w1.write(cell+'\t'+li[2]+'\t'+li[3]+'\n')
    if cluster.startswith('nsclc'):
        for cell in dict_nsclc[cluster]:
            w2.write(cell+'\t'+li[2]+'\t'+li[3]+'\n')
    if cluster.startswith('normal'):
        for cell in dict_normal[cluster]:
            w3.write(cell+'\t'+li[2]+'\t'+li[3]+'\n')

w1.close()
w2.close()
w3.close()

###############calculate mutation rate

def output_count(file_in,file_out):
    a1=open(file_in,'r')
    w=open(file_out,'w')
    dict_cell2count={}
    for line in a1:
        li=line.split()
        if li[0] not in dict_cell2count:
            count=0
            if li[1] !='2':
                count=count+1
                dict_cell2count[li[0]]=count
        else:
            if li[1] != '2':
                dict_cell2count[li[0]]=dict_cell2count[li[0]]+1
    for cell in dict_cell2count:
        w.write(cell+'\t'+str(dict_cell2count[cell])+'\n')
    w.close()

output_count('cells_cnv_normal.txt','cells_cnv_percentage_normal.txt')
output_count('cells_cnv_nsclc.txt','cells_cnv_percentage_nsclc.txt')
output_count('cells_cnv_sclc.txt','cells_cnv_percentage_sclc.txt')

############### draw hist plot
sclc_percentage <- read.table('cells_cnv_percentage_sclc.txt')
nsclc_percentage <- read.table('cells_cnv_percentage_nsclc.txt')
normal_percentage <- read.table('cells_cnv_percentage_normal.txt')

hist(sclc_percentage$V3,  col = rgb(0, 0, 1, 0.3), breaks = 30, freq = F, ylim = c(0,15), border = NA)
hist(nsclc_percentage$V3,  col = rgb(1, 0, 0, 0.3), breaks = 10, freq = F, add = TRUE, border = NA)
hist(normal_percentage$V3,  col = NA, breaks = 10, freq = F, add = TRUE)

abline(v=0.15, col = 'red')

legend("topright", legend = c("NSCLC", "SCLC","Normal"), fill = c("pink", "#B2B0EC","white"))

sclc_cells <- sclc_percentage$V1[sclc_percentage$V3 > 0.15]
write(sclc_cells, file = "sclc_cells.txt")
nsclc_cells <- nsclc_percentage$V1[nsclc_percentage$V3 > 0.15]
write(nsclc_cells, file = "nsclc_cells.txt")

non_tumor_cells <- names(deviations)[!names(deviations) %in% cancer_cells | grepl("^HTA8_3", names(deviations))]
write(non_tumor_cells , file = "non_cancer_cells.txt")


############################################################
        extract the names of cells in ref type
############################################################
setwd('/home/yq238/project/scRNA_seq/HTAN2/09_epithelial_subtype')
ref <- readRDS('/home/yq238/project/scRNA_seq/HTAN2/09_epithelial_subtype/ref_epi.rds')

write.table(rownames(ref@meta.data[ref@meta.data$author_cell_type=='NSCLC',]),'NSCLC_cell_ref.txt',sep = '\t',quote = F, col.names = F, row.names = F)

selected_cells <- rownames(ref@meta.data[!ref@meta.data$author_cell_type %in% c("SCLC-P", "SCLC-A", "SCLC-N", "NSCLC"),])
write.table(selected_cells, 'non_tumor_cells_ref.txt', sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

############################################################
                 add annotation
############################################################
load("/home/yq238/project/scRNA_seq/HTAN2/09_epithelial_subtype_revised/seurat_epithelial_ne_annotation.RData")

tumor_cells <- read.table("cancer_cells.txt")
tumor_cells <- tumor_cells$V1


seurat_epithelial@meta.data$tumorstatus <- ifelse(rownames(seurat_epithelial@meta.data) %in% tumor_cells, "cancer",
 "normal")

save(seurat_epithelial, file="seurat_epithelial_tumorstatus_annotation.RData")


pdf('tumors_cells_annotation_umap.pdf', width = 7 , height = 6)
DimPlot(seurat_epithelial, reduction = "umap.mnn", label = TRUE, repel = TRUE, group.by='tumorstatus')
dev.off()

write.table(rownames(seurat_epithelial@meta.data[seurat_epithelial@meta.data$tumorstatus=='normal',]),'normal_cells.txt',sep = '\t',quote = F, col.names = F, row.names = F)



