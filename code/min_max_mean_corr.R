

####################### mean, min and max correlation  #################################

cor_result <- get(load("../rdas/ash_cor_PD_gtex_tissues.rda"))
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

david_chart <- read.table("../utilities/david_fucntional_chart.txt")
