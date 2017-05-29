

###############  gene wide CorShrink  (tissue pair wise)  ###################

cor_data <- get(load("../output/cor_tissues_non_ash_voom_pearson.rda"))
corvec_data <- apply(cor_data, 3, function(x) return(x[lower.tri(x)]))
common_samples <- get(load("../output/common_samples.rda"))
common_samples_vec <- common_samples[lower.tri(common_samples)]
library(ashr)

ash_out_mat <- matrix(0, dim(corvec_data)[1], dim(corvec_data)[2])
for(j in 1:dim(ash_out_mat)[1]){
  ash_out_mat[j,] <- CorShrinkMLvec(corvec_data[1,], common_samples_vec[1], sd_boot = FALSE,
                                    ash.control = list(mixcompdist = "halfuniform",
                                                       mode = "estimate", nullweight = 100))
}

save(ash_out_mat, file = "../output/genewide_ash_out_mat.rda")

ash_tissue_mat <- array(0, c(dim(cor_data)[1], dim(cor_data)[2], dim(cor_data)[3]))
for(m in 1:dim(tissue_mat)[3]){
  temp_mat <- matrix(0, dim(cor_data)[1], dim(cor_data)[2])
  temp_mat[lower.tri(temp_mat)] <- ash_out_mat[,m]
  temp_mat_2 <- temp_mat + t(temp_mat)
  diag(temp_mat_2) <- 1
  ash_tissue_mat[,,m] <- temp_mat_2
}

save(ash_tissue_mat, file = "../output/genewide_ash_out_tissue_mat.rda")
