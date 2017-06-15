
##########  FLASH on the age data  #################################

data <- get(load("../output/age_cor_nonmash.rda"))
betahat <- data$betahat
sebetahat <- data$sebetahat
zdat <- betahat/sebetahat
out <- greedy(zdat, K=3, flash_para = list(partype = "known", sigmae2_true = rep(1, dim(zdat)[2])))


set.seed(99)
N = 100
P = 200
data = sim_hd(N, P, SF = 0.5, SL=0.5, signal =1, a = rep(1,N), b = rchisq(P,1), mu = 0)
Y = data$Y
E = data$Error
gf_3 = flash(Y,objtype = "l", partype = "known",sigmae2_true = data$sig2_true[1,])
gf_4 = flash(Y,objtype = "l", partype = "known",sigmae2_true = rep(1, dim(Y)[2]))

out <- greedy(Y, K=3, flash_para = list(partype = "known", sigmae2_true = rep(1, dim(Y)[2])))
ba164c770277aef192430f4883fcc18f12b61c3b


library(ashr)
library(flashr)

data <- get(load("../output/sex_cor_nonmash.rda"))
betahat <- data$betahat
sebetahat <- data$sebetahat
zdat <- betahat/sebetahat
out <- greedy(zdat, K=10, flash_para = list(partype = "known", sigmae2_true = rep(1, dim(zdat)[2])))

save(out, file = "../output/sex_cor_flash.rda")

library(ashr)
library(flashr)

data <- get(load("../output/suicide_cor_nonmash.rda"))
betahat <- data$betahat
sebetahat <- data$sebetahat
zdat <- betahat/sebetahat
out <- greedy(zdat, K=10, flash_para = list(partype = "known", sigmae2_true = rep(1, dim(zdat)[2])))

save(out, file = "../output/suicide_cor_flash.rda")



