
#########################   flash systemic post processing   ########################

flash_out <- get(load("../output/sex_cor_flash.rda"))

common_samples <- get(load("../output/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

#library(data.table)
#data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"))
#matdata <- t(data[,-c(1,2)])

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))


U <- unique(tissue_labels)

annotation <- data.frame(sample_id = 1:NROW(flash_out$l),
                         label = U)
rownames(flash_out$l) <- 1:dim(flash_out$l)[1]
FactorGGBar(flash_out$l, annotation = annotation)



flash_out <- get(load("../output/circadian_cor_flash.rda"))
annotation <- data.frame(sample_id = 1:NROW(flash_out$l),
                         label = U)
rownames(flash_out$l) <- 1:dim(flash_out$l)[1]
FactorGGBar(flash_out$l, annotation = annotation)


flash_out <- get(load("../output/bmi_cor_flash.rda"))
annotation <- data.frame(sample_id = 1:NROW(flash_out$l),
                         label = U)
rownames(flash_out$l) <- 1:dim(flash_out$l)[1]
FactorGGBar(flash_out$l, annotation = annotation)

flash_out <- get(load("../output/age_cor_flash.rda"))
annotation <- data.frame(sample_id = 1:NROW(flash_out$l),
                         label = U)
rownames(flash_out$l) <- 1:dim(flash_out$l)[1]
FactorGGBar(flash_out$l, annotation = annotation)

flash_out <- get(load("../output/suicide_cor_flash.rda"))
annotation <- data.frame(sample_id = 1:NROW(flash_out$l),
                         label = U)
rownames(flash_out$l) <- 1:dim(flash_out$l)[1]
FactorGGBar(flash_out$l, annotation = annotation)
