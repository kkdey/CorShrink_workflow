

####################   cellcycler  prepare data  #########################

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

liver_data  <- person_tissue_genes_1[,which(U=="Liver"),]


###################  liver specific genes  ###########################

liver_specific_genes <- read.table("../utilities/circadian_mash_3/liver.txt")[,1]

liver_data_2 <- liver_data[,match(as.character(liver_specific_genes), gene_names_1)]
colnames(liver_data_2) <- as.character(liver_specific_genes)


plot(uncos_times_1, liver_data_2[,200], pch=20, cex=1)

liver_data_3 <- liver_data_2[which(!is.na(liver_data_2[,1])),]
dim(liver_data_3)

uncos_times_2 <- uncos_times_1[which(!is.na(liver_data_2[,1]))]

mash_dat <- list("data" = liver_data_3,
                 "times" = uncos_times_2)
save(mash_dat, file = "../output/cellcyler_data_liver_mash.rda")

index <- c(1, 10, 20)
plot(uncos_times_2, liver_data_3[,index[1]])
plot(uncos_times_2, liver_data_3[,index[2]])
plot(uncos_times_2, liver_data_3[,index[3]])

plot(uncos_times_2, liver_data_3[,200])


indices <- match("ENSG00000140057", gene_names_1)
par(mfrow=c(3,3))
for(j in 1:53){
  plot(uncos_times_1, person_tissue_genes_1[,j, indices],
       ylab = dimnames(person_tissue_genes)[[2]][j])
}


######################   circadian genes   ##############################

pathways <- read.delim("../utilities/CPDB_pathways_genes.tab")
pathway_names <- pathways[,1]
grep("Circadian",pathway_names)
circadian_pathways <- pathways[grep("Circadian", pathway_names),]
circadian_genes <- list()
for (i in 1:dim(circadian_pathways)[1]){
  circadian_genes[[i]] <- strsplit(as.character(circadian_pathways[i,4]), "[,]")[[1]]
}

circ_names <- Reduce(union, circadian_genes)

indices1 <- match(as.character(circ_names), gene_names_1)
indices2 <- indices1[!is.na(indices1)]

liver_data_circ <- liver_data[,indices2]
colnames(liver_data_circ) <- as.character(circ_names[!is.na(indices1)])

liver_data_circ_2 <- liver_data_circ[which(!is.na(liver_data_circ[,1])),]
dim(liver_data_circ_2)

uncos_times_circ <- uncos_times_1[which(!is.na(liver_data_circ[,1]))]

circ_dat <- list("data" = liver_data_circ_2,
                 "times" = uncos_times_circ)
save(circ_dat, file = "../output/cellcyler_data_liver_circadian.rda")

plot(uncos_times_circ, liver_data_circ_2[,1])
plot(uncos_times_circ, liver_data_circ_2[,10])
plot(uncos_times_circ, liver_data_circ_2[,50])
