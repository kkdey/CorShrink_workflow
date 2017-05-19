
######################      corpcor vs corshrink     ############################

person_label=read.table("../data/GTEX_V6/person_identifier_labels_with_numbers.txt");
samples_id <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,1]

samples_person <- sapply(samples_id, function(x) return(paste0(strsplit(as.character(x), "-")[[1]][1:2], collapse ="-")))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

xtabs(~ samples_person + tissue_labels)

unique_persons <- unique(samples_person)
unique_tissues <- unique(tissue_labels)

library(data.table)
data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"))
matdata <- (data[,-c(1,2)])
voom_matdata <- t(limma::voom(matdata)$E)

mat <- matrix(0, length(unique_persons), length(unique_tissues))


person_tissue_genes <- array(NA, c(length(unique_persons), length(unique_tissues), dim(voom_matdata)[2]))

for(u in 1:length(unique_persons)){
    index_samp <- which(samples_person == unique_persons[u])
    tissue_lab_samp <- tissue_labels[index_samp]
    person_tissue_genes[u, match(tissue_lab_samp, unique_tissues), ] = voom_matdata[index_samp, ]
    cat("We are at person : ", u, "\n")
}

save(person_tissue_genes, file = "../rdas/person_tissue_genes_voom.rda")



#############  exmaple corpcor vs corshrink for a single gene  ############################

numg <- 100
mat <- person_tissue_genes[,,numg]
rank_mat <- matrix(0, dim(mat)[1], dim(mat)[2])
for(u in 1:dim(mat)[1]){
  temp <- mat[u, ]
  w <- which(temp == 0);
  mat[u, w] = sample(temp[temp != 0], length(w), replace=TRUE)
  rank_mat[u,] <- rank(mat[u,])
}



cov_mat <- cov(rank_mat);
system.time(strimmer_sample <- corpcor::cov.shrink(mat))
system.time(glasso_sample_005 <- glasso::glasso(cov_mat, rho = 0.05))
system.time(glasso_sample_05 <- glasso::glasso(cov_mat, rho = 0.5))
system.time(glasso_sample_1 <- glasso::glasso(cov_mat, rho = 1))
system.time(glasso_sample_10 <- glasso::glasso(cov_mat, rho = 10))
system.time(cov_sample_ML <-  CorShrink::CorShrinkML(cov2cor(cov_mat), nsamp_mat = 550, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal",
                                                             nullweight = 1)))


order_index <- c();
U <- unique_tissues
order_index <- c(order_index, grep("Brain", U))
order_index <- c(order_index, grep("fibroblast", U))
order_index <- c(order_index, grep("EBV", U))
order_index <- c(order_index, grep("Spleen", U))
order_index <- c(order_index, grep("Whole Blood", U))
order_index <- c(order_index, grep("Muscle - Skeletal", U))
order_index <- c(order_index, grep("Liver", U))
order_index <- c(order_index, grep("Pancreas", U))
order_index <- c(order_index, grep("Stomach", U))
order_index <- c(order_index, grep("Kidney - Cortex", U))
order_index <- c(order_index, grep("Adrenal Gland", U))
order_index <- c(order_index, grep("Colon - Transverse", U))
order_index <- c(order_index, grep("Small Intestine - Terminal Ileum", U))
order_index <- c(order_index, grep("Heart - Atrial Appendage", U))
order_index <- c(order_index, grep("Heart - Left Ventricle", U))
order_index <- c(order_index, grep("Minor Salivary Gland", U))
order_index <- c(order_index, grep("Skin - Sun Exposed", U))
order_index <- c(order_index, grep("Skin - Not Sun Exposed", U))
order_index <- c(order_index, grep("Lung", U))
order_index <- c(order_index, grep("Ovary", U))
order_index <- c(order_index, grep("Thyroid", U))
order_index <- c(order_index, grep("Pituitary", U))
order_index <- c(order_index, grep("Testis", U))
order_index <- c(order_index, grep("Nerve - Tibial", U))
order_index <- c(order_index, grep("Breast - Mammary Tissue", U))
order_index <- c(order_index, grep("Adipose - Visceral", U))
order_index <- c(order_index, grep("Adipose - Subcutaneous", U))
order_index <- c(order_index, grep("Artery - Coronary", U))
order_index <- c(order_index, grep("Artery - Tibial", U))
order_index <- c(order_index, grep("Artery - Aorta", U))
order_index <- c(order_index, grep("Esophagus - Mucosa", U))
order_index <- c(order_index, grep("Vagina", U))
order_index <- c(order_index, grep("Cervix - Endocervix", U))
order_index <- c(order_index, grep("Esophagus - Gastroesophageal Junction", U))
order_index <- c(order_index, grep("Colon - Sigmoid", U))
order_index <- c(order_index, grep("Esophagus - Muscularis", U))
order_index <- c(order_index, grep("Cervix - Ectocervix", U))
order_index <- c(order_index, grep("Fallopian", U))
order_index <- c(order_index, grep("Prostate", U))
order_index <- c(order_index, grep("Uterus", U))
order_index <- c(order_index, grep("Bladder", U))


cor_result <- get(load("../rdas/ash_cor_only_gtex_tissues.rda"))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[order_index,order_index, numg]),
      col=col, main=paste0("corshrink: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)
axis(2, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cov2cor(strimmer_sample))[order_index, order_index],
      col=col, main=paste0("shafer strimmer: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)
axis(2, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cov2cor(glasso_sample_005$w))[order_index, order_index],
      col=col, main=paste0("glasso 0.05: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)
axis(2, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cov2cor(glasso_sample_05$w))[order_index, order_index],
      col=col, main=paste0("glasso 0.5: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)
axis(2, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cov_sample_ML$ash_cor_only)[order_index, order_index],
      col=col, main=paste0("corshrink on fake: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)
axis(2, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)
