

##############  genewide ash on the suicide data  ####################

library(ashr)

person_tissue_genes <- get(load("../output/person_tissue_genes_voom.rda"))

data <- get(load("../output/cva_v_anoxia_cor_nonmash.rda"))
betahat <- data$betahat
sebetahat <- data$sebetahat

pos_mean <- matrix(0, dim(betahat)[1], dim(betahat)[2])
pos_sd <- matrix(0, dim(betahat)[1], dim(betahat)[2])
pos_lfsr <- matrix(0, dim(betahat)[1], dim(betahat)[2])

for(l in 1:dim(betahat)[2]){
  ash_fit <- ashr::ash(betahat[,l], sebetahat[,l], mixcompdist = "halfuniform",
                       nullweight = 100)
  pos_mean[,l] <- ash_fit$result$PosteriorMean
  pos_sd[,l] <- ash_fit$result$PosteriorSD
  pos_lfsr[,l] <- ash_fit$result$lfsr
  cat("We are at iteration : ", l, "\n")
}

rownames(pos_mean) <- dimnames(person_tissue_genes)[[2]]
rownames(pos_sd) <- dimnames(person_tissue_genes)[[2]]
rownames(pos_lfsr) <- dimnames(person_tissue_genes)[[2]]


colnames(pos_mean) <- dimnames(person_tissue_genes)[[3]]
colnames(pos_sd) <- dimnames(person_tissue_genes)[[3]]
colnames(pos_lfsr) <- dimnames(person_tissue_genes)[[3]]



dim(pos_lfsr)

ll <- list("post.mean" = pos_mean,
           "post.sd" = pos_sd,
           "post.lfsr" = pos_lfsr)

save(ll, file = "../output/cva_v_anoxia_cor_ash.rda")

plot(betahat[10,], col="red")
points(ll$post.mean[10,], col="blue")

num_lfsr_below_thresh_1 <- apply(ll$post.lfsr[grep("Brain", rownames(ll$post.lfsr)),], 2, function(x) return(length(which(x < 0.1))))
num_lfsr_below_thresh_2 <- apply(ll$post.lfsr[-grep("Brain", rownames(ll$post.lfsr)),], 2, function(x) return(length(which(x < 0.1))))
diff <- num_lfsr_below_thresh_1 - num_lfsr_below_thresh_2

imp_genes_1 <- order(diff, decreasing = TRUE)[1:10]
names <- colnames(pos_mean)[imp_genes_1]


num_lfsr_below_thresh <- apply(ll$post.lfsr, 2, function(x) return(length(which(x < 0.1))))

imp_genes <- order(num_lfsr_below_thresh, decreasing = TRUE)[1:25]
num_lfsr_below_thresh[imp_genes]

names <- colnames(pos_mean)[imp_genes]
write.table(names, file = "../utilities/cva_heart_genes_ash.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

grep("Brain", rownames(ll$post.lfsr))
