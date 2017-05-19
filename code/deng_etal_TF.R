

##################  CorShrink on Deng et al TF data   ################################

deng_tf_data <- get(load("../ShrinkCovpaper/deng_counts_tf.rda"))
non_zero_indices <- apply(deng_tf_data, 1, function(x) length(which(x!=0)))

deng_tf_data <- deng_tf_data[which(non_zero_indices > 30), ]

voom_deng_tf <- limma::voom(deng_tf_data)$E
cov_sample <- cov(t(deng_tf_data))

library(reshape2)
cor_sample_melt <- melt(cov2cor(cov_sample))
mod_cor_sample_melt <- cor_sample_melt[cor_sample_melt[,1] != cor_sample_melt[,2],]

mod_cor_sample_melt[which.max(mod_cor_sample_melt[,3]),1]
mod_cor_sample_melt[which.max(mod_cor_sample_melt[,3]),2]

library(devtools)

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()

deng.meta_data <- Biobase::pData(Deng2014MouseESC)

cell_meta <- deng.meta_data$cell_type[which(non_zero_indices > 30)]
index <- grep(paste(mod_cor_sample_melt[which.max(mod_cor_sample_melt[,3]),2]), rownames(deng_tf_data))
plot(1:dim(deng_tf_data)[2], log(deng_tf_data[index,]+1), type="l", col="red",  ylab="log expr.", main=paste(mod_cor_sample_melt[which.max(mod_cor_sample_melt[,3]),2],"- expression"), xaxt="n", xlab="")

labels = match(unique(cell_meta), cell_meta);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(deng_tf_data)[2]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(cell_meta),las=3, cex.axis=0.8);


mat <- cov2cor(cov_sample)
x <- as.matrix(mat)
library(gplots)
png(file = "heatmap.png")
heatmap.2(x, trace="none", symkey=TRUE, key=TRUE,
          labRow = FALSE, labCol = FALSE)
dev.off()


system.time(covshrink_ML <- CorShrink::CovShrink((voom_deng_tf), 250, type="ML"));

system.time(corpcor::cov.shrink((voom_deng_tf)));


save(covshrink_ML, file="../ShrinkCovpaper/covshrink_deng_tf.rda")

covshrink_ML <- get(load("../ShrinkCovpaper/covshrink_deng_tf.rda"))

mat <- cov2cor(as.matrix(covshrink_ML))
x <- as.matrix(mat)
library(gplots)
png(file = "heatmap2.png")
heatmap.2(x, trace="none", key=TRUE,
          labRow = FALSE, labCol = FALSE)
dev.off()


covshrink_ML <- get(load("../ShrinkCovpaper/covshrink_deng_tf.rda"))
corshrink_ML <- cov2cor(as.matrix(covshrink_ML))

indices <- which(corshrink_ML > 0.9 & corshrink_ML < 1, arr.ind = TRUE)
rownames(deng_tf_data)[unique(indices[,2])]

write.table(rownames(deng_tf_data)[unique(indices[,2])], "../high_cor_genes_deng_tf.txt", col.names = FALSE, row.names=FALSE, quote=FALSE);
write.table(rownames(deng_tf_data), "../all_genes_deng_tf.txt", col.names = FALSE, row.names=FALSE, quote=FALSE);

num <- 4
index <- grep(paste(rownames(deng_tf_data)[unique(indices[,2])][num]), rownames(deng_tf_data))
plot(1:dim(deng_tf_data)[2], log(deng_tf_data[index,]+1), type="l", col="red",  ylab="log expr.", main=paste(rownames(deng_tf_data)[unique(indices[,2])][num],"- expression"), xaxt="n", xlab="")

labels = match(unique(cell_meta), cell_meta);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(deng_tf_data)[2]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(cell_meta),las=3, cex.axis=0.8);



index <- grep("", rownames(reads))[2]
plot(1:dim(reads)[2], log(reads[index,]+1), type="l", col="red",  ylab="log expr.", main="Tcl1 expression", xaxt="n", xlab="")

labels = match(unique(cell_meta), cell_meta);
abline(v=labels)

labels_low=labels;
labels_up=c(labels[2:length(labels)],dim(reads)[2]);
mid_point=labels_low +0.5*(labels_up-labels_low);

axis(1,at=mid_point, unique(cell_meta),las=3, cex.axis=0.8);

cluster     = hclust(as.dist(1- cov2cor(as.matrix(covshrink_ML))))
dendrogram  = as.dendrogram(cluster)
order.dendrogram(dendrogram)
dendrogram


new_indices <- tail(order.dendrogram(dendrogram), 100);



cor_sample <- cov2cor(cov_sample)
indices <- which(cor_sample > 0.9 & cor_sample < 1, arr.ind = TRUE)
rownames(deng_tf_data)[unique(indices[,2])]

heatmap(cov_sample)

