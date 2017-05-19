

##############  GSEA  corshrink brain  (CountClust)  ######################


cor_result <- get(load("../rdas/ash_cor_only_gtex_tissues.rda"))
common_samples <- get(load("../rdas/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

#library(data.table)
#data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"))
#matdata <- t(data[,-c(1,2)])

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))


U <- unique(tissue_labels)


order_index <- c();
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

tab <- array(0, dim(cor_result)[3])
for(m in 1:dim(cor_result)[3]){
  z <- as.matrix(cor_result[order_index[1:13],order_index[1:13],m])
  vec_z <- z[row(z) > col(z)]
  tab[m] <- length(which(vec_z > 0.2))
}

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[order_index[1:13],order_index[1:13],which.max(tab)]),
      col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
axis(2, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[order_index[1:13],order_index[1:13],which.min(tab)]),
      col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
axis(2, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)


tval <- (tab - mean(tab))/sd(tab)
names(tval) <- gene_names_1

cluster_list <- vector(mode = "list", length = 6)
for(l in 1:6){
  cluster_list[[l]] <- as.character(read.table(paste0("../utilities/gene_names_brain_clus_", l, ".txt"))[,1])
}

library(fgsea)
out <- fgsea(pathways = cluster_list,
             stats = tval,
             nperm = 100000)

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[order_index[1:13],order_index[1:13], grep("ENSG00000059758", gene_names_1)]),
      col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
axis(2, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)


tab <- array(0, dim(cor_result)[3])
for(m in 1:dim(cor_result)[3]){
  z <- as.matrix(cor_result[order_index[1:13],order_index[1:13],m])
  vec_z <- z[row(z) > col(z)]
  tab[m] <- quantile(vec_z, 0.4)
}

tval <- (tab - mean(tab))/sd(tab)
names(tval) <- gene_names_1

out <- fgsea(pathways = cluster_list,
             stats = tval,
             nperm = 100000)


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[order_index[1:13],order_index[1:13], grep("ENSG00000050165", gene_names_1)]),
      col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
axis(2, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)

k <- 6
genes <- as.character(read.table(paste0("../utilities/gene_names_brain_clus_", k, ".txt"))[1:20,1])
par(mfrow=c(2,2))

for(l in 1:length(genes)){
  col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
        rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
  image(as.matrix(cor_result[order_index[1:13],order_index[1:13],grep(paste0(genes[l]),gene_names_1)]),
        col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
        xaxt = "n", yaxt = "n", zlim=c(-1,1))
  axis(1, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
  axis(2, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
}

brain_related_genes <- c();
for(l in 1:length(cluster_list)){
  brain_related_genes <- c(brain_related_genes, cluster_list[[l]]);
}
brain_related_genes <- unique(brain_related_genes)

other_genes <- setdiff(gene_names_1, brain_related_genes)

cluster_list <- vector(mode = "list", length = 6)
for(l in 1:6){
  cluster_list[[l]] <- as.character(read.table(paste0("../utilities/gene_names_brain_clus_", l, ".txt"))[,1])
}

cluster_list[[7]] <- sample(other_genes, 100)

pval_others <- array(0, 50)

for(num in 1:50){
  out <- fgsea(pathways = cluster_list,
               stats = tval,
               nperm = 10000)
  pval_others[num] <- out$padj[7]
  cat("We are at iteration: ", num, "\n")
}

plot(density(pval_others))
