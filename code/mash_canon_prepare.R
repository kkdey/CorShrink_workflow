
library(mashr)
library(ashr)
library(Rmosek)
library(REBayes)

data <- get(load("../output/sex_cor_nonmash.rda"))

betahat <- data$betahat
sebetahat <- data$sebetahat

cat("starting the canonical mash \n")

mash_data <- set_mash_data(t(betahat), t(sebetahat))
U.c = cov_canonical(mash_data)
print(names(U.c))

common_samples <- get(load("../output/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

U <- unique(tissue_labels)

U.adipose1 <- matrix(0, 53, 53)
U.adipose2 <- matrix(0, 53, 53)
l_adipose <- length(grep("Adipose", U))
U.adipose1[grep("Adipose", U), grep("Adipose", U)] <- 0.5*diag(1, l_adipose) + 0.5*rep(1, l_adipose) %*% t(rep(1, l_adipose))
U.adipose2[grep("Adipose", U), grep("Adipose", U)] <- rep(1, l_adipose) %*% t(rep(1, l_adipose))


U.artery1 <- matrix(0, 53, 53)
U.artery2 <- matrix(0, 53, 53)
l_artery <- length(grep("Artery", U))
U.artery1[grep("Artery", U), grep("Artery", U)] <- 0.5*diag(1, l_artery) + 0.5*rep(1, l_artery) %*% t(rep(1, l_artery))
U.artery2[grep("Artery", U), grep("Artery", U)] <- rep(1, l_artery) %*% t(rep(1, l_artery))


U.brain1 <- matrix(0, 53, 53)
U.brain2 <- matrix(0, 53, 53)
l_brain <- length(grep("Brain", U))
U.brain1[grep("Brain", U), grep("Brain", U)] <- 0.5*diag(1, l_brain) + 0.5*rep(1, l_brain) %*% t(rep(1, l_brain))
U.brain2[grep("Brain", U), grep("Brain", U)] <- rep(1, l_brain) %*% t(rep(1, l_brain))


U.nonbrain1 <- matrix(0, 53, 53)
U.nonbrain2 <- matrix(0, 53, 53)
l_nonbrain <- 53 - length(grep("Brain", U))
U.nonbrain1[-grep("Brain", U), -grep("Brain", U)] <- 0.5*diag(1, l_nonbrain) + 0.5*rep(1, l_nonbrain) %*% t(rep(1, l_nonbrain))
U.nonbrain2[-grep("Brain", U), -grep("Brain", U)] <- rep(1, l_nonbrain) %*% t(rep(1, l_nonbrain))


U.cerebellum1 <- matrix(0, 53, 53)
U.cerebellum2 <- matrix(0, 53, 53)
l_cerebellum <- length(grep("Cerebell", U))
U.cerebellum1[grep("Cerebell", U), grep("Cerebell", U)] <- 0.5*diag(1, l_cerebellum) + 0.5*rep(1, l_cerebellum) %*% t(rep(1, l_cerebellum))
U.cerebellum2[grep("Cerebell", U), grep("Cerebell", U)] <- rep(1, l_cerebellum) %*% t(rep(1, l_cerebellum))


U.heart1 <- matrix(0, 53, 53)
U.heart2 <- matrix(0, 53, 53)
l_heart <- length(grep("Heart", U))
U.heart1[grep("Heart", U), grep("Heart", U)] <- 0.5*diag(1, l_heart) + 0.5*rep(1, l_heart) %*% t(rep(1, l_heart))
U.heart2[grep("Heart", U), grep("Heart", U)] <- rep(1, l_heart) %*% t(rep(1, l_heart))


U.skin1 <- matrix(0, 53, 53)
U.skin2 <- matrix(0, 53, 53)
l_skin <- length(grep("Skin", U))
U.skin1[grep("Skin", U), grep("Skin", U)] <- 0.5*diag(1, l_skin) + 0.5*rep(1, l_skin) %*% t(rep(1, l_skin))
U.skin2[grep("Skin", U), grep("Skin", U)] <- rep(1, l_skin) %*% t(rep(1, l_skin))



U.pool <- U.c
U.pool[["adipose_0_5"]] <- U.adipose1
U.pool[["adipose_1"]] <- U.adipose2
U.pool[["artery_0_5"]] <- U.artery1
U.pool[["artery_1"]] <- U.artery2
U.pool[["brain_0_5"]] <- U.brain1
U.pool[["brain_1"]] <- U.brain2
U.pool[["cerebellum_0_5"]] <- U.cerebellum1
U.pool[["cerebellum_1"]] <- U.cerebellum2
U.pool[["heart_0_5"]] <- U.heart1
U.pool[["heart_1"]] <- U.heart2
U.pool[["skin_0_5"]] <- U.skin1
U.pool[["skin_1"]] <- U.skin2
U.pool[["non_brain_0_5"]] <- U.nonbrain1
U.pool[["non_brain_1"]] <- U.nonbrain2


#############  canonical mash performance ###############

m.pool = mash(mash_data, U.pool, optmethod = "mixEM")
save(m.c, file = "../output/mash_sex_canonical_pool_2.rda")


