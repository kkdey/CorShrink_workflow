

##############  GTEx all tissues ash correlation  ###########################

library(ashr)
library(CorShrink)

cor_result <- get(load("../rdas/cor_tissues_non_ash.rda"))
common_samples <- get(load("../rdas/common_samples.rda"))

ash_cor_PD <- array(0, c(dim(cor_result)[1], dim(cor_result)[2], dim(cor_result)[3]))
ash_cor_only <- array(0, c(dim(cor_result)[1], dim(cor_result)[2], dim(cor_result)[3]))

for(num in 1:2){
  cor_mat <- diag(1,53)+cor_result[,,num]
  system.time(cor_sample_ML <-  CorShrinkML(cor_mat, common_samples, sd_boot = FALSE,
                                            ash.control = list(mixcompdist = "normal",
                                                               nullweight = 100)))
  ash_cor_PD[,,num] <- as.matrix(cor_sample_ML$ash_cor_PD)
  ash_cor_only[,,num] <- as.matrix(cor_sample_ML$ash_cor_only)
}

save(ash_cor_only, file = "../rdas/ash_cor_only_gtex_tissues.rda")
save(ash_cor_PD, file = "../rdas/ash_cor_PD_gtex_tissues.rda")




