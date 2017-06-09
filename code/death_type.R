

############## Death type analysis  ######################

sample_attributes <- get(load("../output/sample_attributes_filtered.rda"))

manner_of_death <- sample_attributes$DTHMNNR
cause_of_death <- sample_attributes$DTHCOD

#suicide_indices <- grep("Suicide", manner_of_death)
accident_indices <- grep("Accident", manner_of_death)
natural_indices <- grep("Natural", manner_of_death)

total_indices <- c(accident_indices, natural_indices)

person_tissue_genes <- get(load("../output/person_tissue_genes_voom.rda"))

person_tissue_genes_1 <- person_tissue_genes[total_indices,,]

fac_manner_death <- c(rep("Accident", length(accident_indices)),
                      rep("Natural", length(natural_indices)))

betahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
sebetahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
pval_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
rsquare_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])



for(i in 1:dim(person_tissue_genes)[2]){
  for(j in 1:dim(person_tissue_genes)[3]){
    x <- person_tissue_genes_1[, i, j]
    expr <- x[!is.na(x)]
    if(length(expr) <= 3){
      betahat <- 0
      sebetahat <- 1
      pval <- 0
      rsquare <- 0
    }else if (length(unique(fac_manner_death[!is.na(x)])) < 2){
      betahat <- 0
      sebetahat <- 1
      pval <- 0
      rsquare <- 0
    }else{
      dtype <- fac_manner_death[!is.na(x)]
      fit <- lm(expr ~ factor(dtype))
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
save(ll, file = "../output/accident_cor_nonmash.rda")



