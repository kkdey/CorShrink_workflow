

#################  prepare the FLASH data  ############################

data <- get(load("../rdas/ash_cor_only_voom_pearson_gtex_tissues.rda"))

flash_data <- matrix(0, dim(data)[1]*(dim(data)[1]-1)/2, dim(data)[3])

for(g in 1:dim(flash_data)[2]){
  z <- data[,,g]
  flash_data[,g] <- z[lower.tri(z)]
}

zbar <- matrix(0, dim(data)[1], dim(data)[2])
zbar[lower.tri(zbar)] <- z[lower.tri(z)]

flash_data_trans <- 0.5*log((1+flash_data)/(1-flash_data))


####################  SVD of the flash data  ###########################

svd_out <- prcomp(flash_data_trans)

u_mat <- svd_out$x
corr_mat_vectors_facs <- u_mat[,1:20]
corr_mat_vectors_facs_trans <- (exp(2*corr_mat_vectors_facs)-1)/(exp(2*corr_mat_vectors_facs)+1)


corr_mat_fac <- array(0, c(dim(data)[1], dim(data)[2], 20));

for(k in 1:20){
  mat2 <- matrix(0, dim(data)[1], dim(data)[2])
  mat2[lower.tri(mat2)] <- corr_mat_vectors_facs_trans[,k]
  diag(mat2) <- 1
  corr_mat_fac[,,k] <- (mat2 + t(mat2))
}


num <- 20
col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))
image(as.matrix(corr_mat_fac[order_index,order_index, num]),
      col=col, main=paste0("sample corr with NA: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
axis(1, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)
axis(2, at = seq(0, 1, length.out = 53), labels = unique_tissues[order_index], las=2, cex.axis = 1)




set.seed(99)
N = 100
P = 600
data = sim_K(K=10,N, P, SF = 0.9, SL = 0.8, signal = 1,noise = 1)
Y = data$Y
E = data$Error
plot(svd(Y)$d)
ggd = greedy(Y,K = 100)
#eigen_plot(ggd,sum(Y^2))
gbf <- backfitting(Y, ggd, maxiter_bf = 10, parallel = TRUE)
