

###############  Deng et al 2014 analysis   ##########################


deng_tf_data <- get(load("../rda/deng_counts_tf.rda"))
non_zero_indices <- apply(deng_tf_data, 1, function(x) length(which(x!=0)))

deng_tf_data <- deng_tf_data[which(non_zero_indices > 30), ]
voom_deng_tf <- limma::voom(deng_tf_data)$E

system.time(covshrink_ML <- CorShrink::CovShrink((voom_deng_tf), nsamples = 50,
                                                 sd_boot = FALSE,
                                                 partial = FALSE,
                                                 type="ML"));

save(covshrink_ML, file="../rda/covshrink_deng_tf-noboot.rda")

system.time(covshrink_ML_boot <- CorShrink::CovShrink((voom_deng_tf), nsamples = 50,
                                                 sd_boot = TRUE,
                                                 partial = FALSE,
                                                 type="ML"));

save(covshrink_ML_boot, file="../rda/covshrink_deng_tf-boot.rda")


library(fields)
set.seed(1)
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

cov_sample <- cov(voom_deng_tf)

num <- 100
cols = colorRampPalette(c("red", "white", "blue"))(100)
image.plot(cov2cor(cov_sample), col=cols, nlevel=num, main="sample.corr", cex.main=1)
image.plot(cov2cor(as.matrix(covshrink_ML)), col=cols, nlevel=num, main="CorShrink.ML", cex.main=1)
