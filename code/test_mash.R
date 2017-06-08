

###########  test mash  #####################

library(mashr)
library(ashr)

data <- get(load("../output/death_time_cor_nonmash.rda"))

betahat <- round(data$betahat,4)
sebetahat <- round(data$sebetahat,4)

cat("starting the canonical mash \n")

mash_data <- set_mash_data(t(betahat), t(sebetahat))
U.c = cov_canonical(mash_data)
print(names(U.c))

#############  canonical mash performance ###############

m.c = mash(mash_data, U.c)
save(m.c, file = "../output/mash_circadian_canonical.rda")

############  strong signals detection #############

m.1by1 = mash_1by1(mash_data)
strong = get_significant_results(m.1by1,0.05)

###########  PCA  Us  #################
U.pca = cov_pca(mash_data,5,strong)
print(names(U.pca))

#############  extreme deconvolution #################

cat("starting the extreme deconvolution mash \n")

U.ed = cov_ed(mash_data, U.pca, strong)
m.ed = mash(mash_data, U.ed)
save(m.ed, file = "../output/mash_circadian_exdeconv.rda")

###########   pooled mash   #######################

cat("starting the pooled mash \n")
U.pool <- c(U.c,U.ed)
m.pool = mash(mash_data, U.pool)
save(m.pool, file = "../output/mash_circadian_pool.rda")




