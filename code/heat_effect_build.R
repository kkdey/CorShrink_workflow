

#################  Heat effect relationship  ####################################

sample_attributes <- get(load("../output/sample_attributes_filtered.rda"))
person_tissue_genes <- get(load("../output/person_tissue_genes_voom.rda"))

betahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
sebetahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
pval_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
rsquare_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])

temp_val <- sample_attributes$TRCRTMP
temp_val_unit <- sample_attributes$TRCRTMPU

temp_val_unit[329] <- "F"

which_celsius <- which(temp_val_unit == "C")

new_temp <- array(0, length(temp_val))
new_temp[which_celsius] <- (temp_val[which_celsius]*9)/5 + 32
new_temp[-which_celsius] <- temp_val[-which_celsius]

na_indices <- union(which(is.na(new_temp)), which(new_temp == 0))

person_tissue_genes_1 <- person_tissue_genes[-na_indices,,]
new_temp_1 <- new_temp[-na_indices]

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
      heat <- new_temp_1[!is.na(x)]
      fit <- lm(expr ~ log(heat))
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

save(ll, file = "../output/heat_cor_nonmash.rda")

