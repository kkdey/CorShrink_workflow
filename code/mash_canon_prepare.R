
library(mashr)
library(ashr)
library(Rmosek)
library(REBayes)

data <- get(load("../output/death_time_cor_nonmash.rda"))

betahat <- data$betahat
sebetahat <- data$sebetahat

cat("starting the canonical mash \n")

mash_data <- set_mash_data(t(betahat), t(sebetahat))
U.c = cov_canonical(mash_data)
print(names(U.c))

common_samples <- get(load("../output/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

U <- unique(tissue_labels)

U.adipose <- matrix(0, 53, 53)
U.adipose[grep("Adipose", U), grep("Adipose", U)] <- 1

U.artery <- matrix(0, 53, 53)
U.artery[grep("Artery", U), grep("Artery", U)] <- 1

U.brain <- matrix(0, 53, 53)
U.brain[grep("Brain", U), grep("Brain", U)] <- 1

U.cerebellum <- matrix(0, 53, 53)
U.cerebellum[grep("Cerebell", U), grep("Cerebell", U)] <- 1

U.heart <- matrix(0, 53, 53)
U.heart[grep("Heart", U), grep("Heart", U)] <- 1

U.skin <- matrix(0, 53, 53)
U.skin[grep("Skin", U), grep("Skin", U)] <- 1


U.pool <- U.c
U.pool[["adipose"]] <- U.adipose
U.pool[["artery"]] <- U.artery
U.pool[["brain"]] <- U.brain
U.pool[["cerebellum"]] <- U.cerebellum
U.pool[["heart"]] <- U.heart
U.pool[["skin"]] <- U.skin


#############  canonical mash performance ###############

m.pool = mash(mash_data, U.pool, optmethod = "mixEM")
save(m.c, file = "../output/mash_circadian_canonical_pool.rda")


