

#############  Time of death relationship (Circadian stuff)  #########################

sample_attributes <- get(load("../output/sample_attributes_filtered.rda"))
person_tissue_genes <- get(load("../output/person_tissue_genes_voom.rda"))

betahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
sebetahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
pval_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
rsquare_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])

times <- strsplit(as.character(sample_attributes$DTHTIME), "[:]")
time_stamp <- unlist(lapply(times, function(x){
  y <- as.numeric(x[1])
  z <- as.numeric(x[2])
  w <- y*60+z
  return(w)
}))/(24*60)

cos_times <- cos(2*pi*time_stamp)
na_indices <- which(is.na(cos_times))

person_tissue_genes_1 <- person_tissue_genes[-na_indices,,]
cos_times_1 <- cos_times[-na_indices]

for(i in 1:dim(person_tissue_genes_1)[2]){
  for(j in 1:dim(person_tissue_genes_1)[3]){
    x <- person_tissue_genes_1[, i, j]
    expr <- x[!is.na(x)]
    if(length(expr) <= 3){
      betahat <- 0
      sebetahat <- 1
      pval <- 0
      rsquare <- 0
    }else{
      cos_t <- cos_times_1[!is.na(x)]
      fit <- lm(expr ~ cos_t)
      betahat <- summary(fit)$coefficients[2,1]
      sebetahat <- summary(fit)$coefficients[2,2]
      pval <- summary(fit)$coefficients[2,4]
      rsquare <- summary(fit)$r.squared
    }

    betahat_mat[i,j] <- betahat
    sebetahat_mat[i,j] <- sebetahat
    pval_mat[i,j] <- pval
    rsquare_mat[i,j] <- rsquare

    cat("We are at iteration with ", i, " and ", j, "\n")
  }
}

ll <- list("betahat" = betahat_mat,
           "sebetahat" = sebetahat_mat,
           "pval" = pval_mat,
           "rsquare" = rsquare_mat)

save(ll, file = "../output/death_time_cor_nonmash.rda")

