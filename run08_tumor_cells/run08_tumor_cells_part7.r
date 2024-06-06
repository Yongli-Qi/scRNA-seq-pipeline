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

cancer_cells <- c(sclc_cells, nsclc_cells)
write(cancer_cells, file = "cancer_cells.txt")