
sample_attributes <- get(load("../output/sample_attributes_filtered.rda"))
person_tissue_genes <- get(load("../output/person_tissue_genes_voom.rda"))

betahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])
sebetahat_mat <- matrix(0, dim(person_tissue_genes)[2], dim(person_tissue_genes)[3])

for(i in 1:dim(person_tissue_genes)[2]){
  for(j in 1:dim(person_tissue_genes)[3]){
    x <- person_tissue_genes[, i, j]
    expr <- x[!is.na(x)]
    if(length(expr) <= 3){
      betahat <- 0
      sebetahat <- 1
    }else{
      height <- sample_attributes$BMI[!is.na(x)]
      fit <- lm(expr ~ height)
      betahat <- summary(fit)$coefficients[2,1]
      sebetahat <- summary(fit)$coefficients[2,2]
    }
    betahat_mat[i,j] <- betahat
    sebetahat_mat[i,j] <- sebetahat
    cat("We are at iteration with ", i, " and ", j, "\n")
  }
}
