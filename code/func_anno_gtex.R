

####################### mean, min and max correlation  #################################

cor_result <- get(load("../rdas/ash_cor_only_gtex_tissues.rda"))
common_samples <- get(load("../rdas/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

#library(data.table)
#data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"))
#matdata <- t(data[,-c(1,2)])

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))


U <- unique(tissue_labels)



cor_result_mean <- apply(cor_result, c(1,2), mean)
cor_result_min <- apply(cor_result, c(1,2), min)
cor_result_max <- apply(cor_result, c(1,2), max)

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result_mean),
      col=col, main=paste0("CorShrink mat: mean "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result_min),
      col=col, main=paste0("CorShrink mat: min "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result_max),
      col=col, main=paste0("CorShrink mat: max "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


####################   brain inhibitory genes    #################################

cor_result_which_min <- apply(cor_result, c(1,2), function(x) return(which.min(x)))

ind <- 6224
col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[,,ind]),
      col=col, main=paste0("CorShrink mat: max "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

gene_names[ind]


##  ribosomal protein S19 pseudogene 3
##  ENSG00000228251: lincRNA
##


####################   Brain inhibitory  genes  #################################

indices_brain <- grep("Brain", U)
cor_result_brain <- cor_result[indices_brain, indices_brain, ]

bad_brain_gene_bool <- apply(cor_result_brain, 3, function(x){
      z = x[row(x) > col(x)]
      if(quantile(z, 0.8) < 0.1) { a = 1;
      return (a)} else {
        a = 0;
        return (a)}
})

bad_brain_gene_indices <- which(bad_brain_gene_bool == 1)
bad_brain_gene_names <- gene_names_1[bad_brain_gene_indices]
bad_brain_gene_names_col <- cbind.data.frame(bad_brain_gene_names)
write.table(bad_brain_gene_names_col, file = "../utilities/bad_brain_gene_names.txt",
            row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(cbind.data.frame(gene_names_1), file = "../utilities/gene_names.txt",
            row.names = FALSE, col.names = FALSE, quote=FALSE)

num <- 240
col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[,,bad_brain_gene_indices[num]]),
      col=col, main=paste0("corshrink: ", bad_brain_gene_names[num]), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

###################  bad brain gene annotation  #################################


out <- mygene::queryMany(bad_brain_gene_names,  scopes="ensembl.gene", fields=c("name", "summary"), species="human");
kable(as.data.frame(out))

grep("tyrosin", as.data.frame(out))

as.data.frame(out)$`Summary`



#####################  read DAVID annotations  ##########################

david_chart <- readLines("../utilities/david_fucntional_chart.txt")
grep("ENSG00000198182", bad_brain_gene_names)


#############   CountClust clusters vs CorShrink plots #####################

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


setdiff(U, U[order_index])

k <- 20
genes <- as.character(read.table(paste0("../utilities/gene_names_clus_", k, ".txt"))[1:6,1])
par(mfrow=c(2,3))

for(l in 1:length(genes)){
  col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
        rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
  image(as.matrix(cor_result[order_index,order_index,grep(paste0(genes[l]),gene_names_1)]),
        col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
        xaxt = "n", yaxt = "n", zlim=c(-1,1))
  axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U[order_index], las=2, cex.axis = 1)
  axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U[order_index], las=2, cex.axis = 1)
}


 k <- 6
 genes <- as.character(read.table(paste0("../utilities/gene_names_brain_clus_", k, ".txt"))[1:6,1])
 par(mfrow=c(2,3))

 for(l in 1:length(genes)){
   col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
         rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
   image(as.matrix(cor_result[order_index[1:13],order_index[1:13],grep(paste0(genes[l]),gene_names_1)]),
         col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
         xaxt = "n", yaxt = "n", zlim=c(-1,1))
   axis(1, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
   axis(2, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
 }

l <- 6

tab <- array(0, dim(cor_result)[3])
for(m in 1:dim(cor_result)[3]){
  z <- as.matrix(cor_result[order_index[1:13],order_index[1:13],m])
  vec_z <- z[row(z) > col(z)]
  tab[m] <- quantile(vec_z, 0.25)
}

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_result[order_index[1:13],order_index[1:13],12314]),
      col=col, main=paste0("corshrink: ", genes[l]), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)
axis(2, at = seq(0, 1, length.out = 13), labels = U[order_index[1:13]], las=2, cex.axis = 1.5)


