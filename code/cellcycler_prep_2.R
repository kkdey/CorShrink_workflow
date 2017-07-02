

################  cellcycleR prepare data part 2   #####################

rm(list=ls())
person_tissue_genes <- get(load("../output/person_tissue_genes_voom.rda"))
sample_attributes <- get(load("../output/sample_attributes_filtered.rda"))

times <- strsplit(as.character(sample_attributes$DTHTIME), "[:]")
time_stamp <- unlist(lapply(times, function(x){
  y <- as.numeric(x[1])
  z <- as.numeric(x[2])
  w <- y*60+z
  return(w)
}))/(24*60)

cos_times <- cos(2*pi*time_stamp)
uncos_times <- 2*pi*time_stamp
na_indices <- which(is.na(cos_times))

person_tissue_genes_1 <- person_tissue_genes[-na_indices,,]
cos_times_1 <- cos_times[-na_indices]
uncos_times_1 <- uncos_times[-na_indices]

common_samples <- get(load("../output/common_samples.rda"))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))

U <- unique(tissue_labels)

full_data <- c();
times <- c()
tissue <- c()

bad_tissues <- c("Cervix - Endocervix", "Cervix - Ectocervix",
                 "Fallopian Tube", "Bladder")

match(bad_tissues, U)

liver_specific_genes <- read.table("../utilities/circadian_mash_3/liver.txt")[,1]

for (i in 1:53){
  if(i != 24 && i!=25 && i !=31 && i!=7){
    temp_data  <- person_tissue_genes_1[,i,]
    temp_data_2 <- temp_data[,match(as.character(liver_specific_genes), gene_names_1)]
    colnames(temp_data_2) <- as.character(liver_specific_genes)
    temp_data_3 <- temp_data_2[which(!is.na(temp_data_2[,1])),]
    full_data <- rbind(full_data, temp_data_3)
    times <- c(times, uncos_times_1[which(!is.na(temp_data_2[,1]))])
    tissue <- c(tissue, rep(as.character(U[i]), dim(temp_data_3)[1]))
  }
}

mash_dat <- list("data" = full_data,
                 "times" = times,
                 "tissue" = tissue)

save(mash_dat, file = "../output/cellcyler_data_all_mash.rda")

#########################  circadian genes   ###############################


pathways <- read.delim("../utilities/CPDB_pathways_genes.tab")
pathway_names <- pathways[,1]
grep("Circadian",pathway_names)
circadian_pathways <- pathways[grep("Circadian", pathway_names),]
circadian_genes <- list()
for (i in 1:dim(circadian_pathways)[1]){
  circadian_genes[[i]] <- strsplit(as.character(circadian_pathways[i,4]), "[,]")[[1]]
}

circ_names <- Reduce(union, circadian_genes)

full_data <- c();
times <- c()
tissue <- c()

bad_tissues <- c("Cervix - Endocervix", "Cervix - Ectocervix",
                 "Fallopian Tube", "Bladder")

match(bad_tissues, U)

indices1 <- match(as.character(circ_names), gene_names_1)
indices2 <- indices1[!is.na(indices1)]
circ_names_2 <- gene_names_1[indices2]

for (i in 1:53){
  if(i != 24 && i!=25 && i !=31 && i!=7){
    temp_data  <- person_tissue_genes_1[,i,]
    temp_data_2 <- temp_data[,match(as.character(circ_names_2), gene_names_1)]
    colnames(temp_data_2) <- as.character(circ_names_2)
    temp_data_3 <- temp_data_2[which(!is.na(temp_data_2[,1])),]
    full_data <- rbind(full_data, temp_data_3)
    times <- c(times, uncos_times_1[which(!is.na(temp_data_2[,1]))])
    tissue <- c(tissue, rep(as.character(U[i]), dim(temp_data_3)[1]))
  }
}

circ_dat <- list("data" = full_data,
                 "times" = times,
                 "tissue" = tissue)
save(circ_dat, file = "../output/cellcyler_data_all_circadian.rda")


