
#####################  person genes tissues JSON in R    #############################

person_tissue_genes <- get(load("../rdas/person_tissue_genes_voom.rda"))

ash_corr <- get(load("../rdas/ash_cor_only_gtex_tissues.rda"))

tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]
U <- unique(tissue_labels)


mat <- ash_corr[,,1]
mat_tab_false <- reshape2::melt(mat)

mat <- round(mat, 2)
rownames(mat) <- U
colnames(mat) <- U

mat_tab <- reshape2::melt(mat)


mat_output_1 <- cbind.data.frame(mat_tab_false[,1:2], mat_tab)

x_mat <- matrix(rnorm(dim(mat_output)[1]*10), dim(mat_output)[1], 10);
colnames(x_mat) <- paste0("X-", 1:dim(x_mat)[2])
y_mat <- matrix(rnorm(dim(mat_output)[1]*10), dim(mat_output)[1], 10);
colnames(y_mat) <- paste0("Y-", 1:dim(y_mat)[2])


colnames(mat_output_1) <- c("row", "col", "tissue1", "tissue2", "value")

mat_output <- cbind.data.frame(mat_output_1, x_mat, y_mat)
write.csv(mat_output, file = "../d3/gene1_corr.csv", quote=FALSE, row.names = FALSE)
