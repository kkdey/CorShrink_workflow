

###################   voom correlation comparison  (corshrink) ######################

cor_result <- get(load("../rdas/cor_tissues_non_ash_voom_pearson.rda"))
common_samples <- get(load("../rdas/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))


U <- unique(tissue_labels)

ash_cor_result <- get(load("../rdas/ash_cor_only_voom_pearson_gtex_tissues.rda"))



num <- grep("ENSG00000244734", gene_names_1)
cor_mat <- diag(1, 53) + cor_result[,,num]

system.time(cor_sample_ML <-  CorShrink::CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal")))


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(cor_mat,
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num]), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

ash_cor_mat <- ash_cor_result[,,num]
col=c(rev(rgb(seq(1,0,length=1000),1, seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_only), col=col,
      main=paste0("Corr mat: ", gene_names_1[num]), cex.main=1, xaxt = "n", yaxt = "n",
      zlim= c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)



num <- grep("ENSG00000132639", gene_names_1)
cor_mat <- diag(1,53)+cor_result[,,num]
col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_mat),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num]), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

system.time(cor_sample_ML <-  CorShrink::CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                                     ash.control = list(mixcompdist = "normal")))


col=c(rev(rgb(seq(1,0,length=1000),1, seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_only), col=col,
      main=paste0("Corr mat: ", gene_names_1[num]), cex.main=1, xaxt = "n", yaxt = "n",
      zlim= c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


######################   voom correlation comparison (spearman) ##########################


cor_result <- get(load("../rdas/cor_tissues_non_ash_voom_spearman.rda"))
common_samples <- get(load("../rdas/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))


U <- unique(tissue_labels)


num <- grep("ENSG00000244734", gene_names_1)
cor_mat <- diag(1, 53) + cor_result[,,num]

system.time(cor_sample_ML <-  CorShrink::CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                                     ash.control = list(mixcompdist = "normal")))


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(cor_mat,
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num]), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

ash_cor_mat <- ash_cor_result[,,num]
col=c(rev(rgb(seq(1,0,length=1000),1, seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_only), col=col,
      main=paste0("Corr mat: ", gene_names_1[num]), cex.main=1, xaxt = "n", yaxt = "n",
      zlim= c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)



num <- grep("ENSG00000132639", gene_names_1)
cor_mat <- diag(1,53)+cor_result[,,num]
col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_mat),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num]), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

system.time(cor_sample_ML <-  CorShrink::CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                                     ash.control = list(mixcompdist = "normal")))


col=c(rev(rgb(seq(1,0,length=1000),1, seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_only), col=col,
      main=paste0("Corr mat: ", gene_names_1[num]), cex.main=1, xaxt = "n", yaxt = "n",
      zlim= c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

#####################  Reads expression vs correlation  ################################


person_tissue_genes <- get(load("../rdas/person_tissue_genes_voom.rda"))
num <- grep("ENSG00000132639", gene_names_1)
hbb_mat <- person_tissue_genes[,,num]

summary(hbb_mat)

library(data.table)
data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"))
matdata <- (data[,-c(1,2)])
voom_matdata <- t(limma::voom(matdata)$E)
plot(matdata[,num], col="red", pch = 20)

mean(matdata[matdata[,num] !=0,num])
mat22 <- person_tissue_genes[,,num]
mean(mat22[mat22!=0])

cor_result[33, 43, num]
ash_cor_result[33, 43, num]

plot(person_tissue_genes[,33,num], person_tissue_genes[,43,num])

cor(person_tissue_genes[,33,num], person_tissue_genes[,43,num], use="complete.obs")
