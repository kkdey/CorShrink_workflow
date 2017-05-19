

###################  GTEx  ash corr  ##################################

cor_result <- get(load("../rda/cor_tissues_non_ash.rda"))
common_samples <- get(load("../rda/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

#library(data.table)
#data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"))
#matdata <- t(data[,-c(1,2)])

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))


U <- unique(tissue_labels)
#rownames(common_samples) <- U
#colnames(common_samples) <- U

#save(common_samples, file = "../rda/common_samples.rda")


# 
# 
# par(mar = par("mar") + c(2, 2, 0, 0))
# col=c(rgb(seq(0,1,length=30),1,seq(0,1,length=30)),rgb(1,seq(1,0,length=30),seq(1,0,length=30)),
#       rgb(seq(1,0,length=30),seq(1,0,length=30), 1))
# image(common_samples, col=col, main="common samples", cex.main=1, xaxt = "n", yaxt = "n")
# axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
# axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
# 

num <- 25
cor_mat <- diag(1,53)+cor_result[,,num]
system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                        ash.control = list(mixcompdist = "normal")))


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num]), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

col=c(rev(rgb(seq(1,0,length=1000),1, seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_mat), col=col,
      main=paste0("Corr mat: ", gene_names_1[num]), cex.main=1, xaxt = "n", yaxt = "n",
      zlim= c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

###########  other shrinkage methods  #########################

num <- 50
cor_mat <- diag(1,53)+cor_result[,,num]
cor_mat_2 <- cor_mat
cor_mat_2[is.na(cor_mat_2)] = 0.0001
#cor_mat_2[(cor_mat_2 == 0)] = 0.01

#system.time(glasso_sample_005 <- glasso::glasso(cor_mat_2, rho = 0.05))
system.time(glasso_sample_05 <- glasso::glasso(cor_mat_2, rho = 0.5))
system.time(glasso_sample_1 <- glasso::glasso(cor_mat_2, rho = 1))
system.time(glasso_sample_10 <- glasso::glasso(cor_mat_2, rho = 10))


col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(glasso_sample_10$w,
      col=col, main="ash corr matrix", cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(glasso_sample_1$w,
      col=col, main="ash corr matrix", cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


######################  different types of ash (nullweight)   ######################################

num <- 5
cor_mat <- diag(1,53)+cor_result[,,num]

system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal",
                                                             nullweight = 100)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " null: ", 100), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal",
                                                             nullweight = 1000)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " null: ", 1000), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal",
                                                             nullweight = 10)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " null: ", 10), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal",
                                                             nullweight = 1)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " null: ", 1), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


######################  different types of ash (uniform)  ###############################
num <- 5
cor_mat <- diag(1,53)+cor_result[,,num]

system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "uniform",
                                                             nullweight = 10)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " comp: uniform "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "+uniform",
                                                             nullweight = 10)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col,main=paste0("CorShrink mat: ", gene_names_1[num], " comp: +uniform "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal",
                                                             nullweight = 10)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " comp: normal "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "-uniform",
                                                             nullweight = 10)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " comp: -uniform "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "halfuniform",
                                                             nullweight = 10,
                                                             mode = 0)))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: ", gene_names_1[num], " comp: halfuniform "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)


system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                          ash.control = list(mixcompdist = "normal",
                                                             nullweight = 10,
                                                             mode = "estimate")))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(cor_sample_ML$ash_cor_PD),
      col=col, main=paste0("CorShrink mat: shrink to est "), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)
axis(2, at = seq(0, 1, length.out = ncol(cor_result[,,1])), labels = U, las=2, cex.axis = 0.4)

