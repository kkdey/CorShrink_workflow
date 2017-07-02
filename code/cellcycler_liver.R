

#############  cellcycler on the liver cells  #####################

data <- get(load("../output/cellcyler_data_liver_mash.rda"))

times <- data$times
expr <- data$data

library(cellcycleR)
out <- np_cell_ordering_class(expr, celltime_levels = 100,
                               num_iter=100,
                               save_path="../output/cell_order_liver.rda")
cell_order_full <- cell_ordering_full(out$signal_intensity, dim(expr)[2])

cycle_data <- expr

G <- dim(cycle_data)[2];
numcells <- dim(cycle_data)[1];
celltime_levels = 100

celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
cell_times_init <- times


cell_times_iter <- times
g <- 200
freq <- 0.5
temp <- rescale(cycle_data[,g], to = c(0, 1))
fit <- lm(temp  ~ cos(freq*cell_times_iter) + sin(freq*cell_times_iter) -1);

plot(times, temp)
points(times, fit$fitted.values, col="blue")
sum(fit$residuals^2)













plot(times, out$cell_times, pch=20, col="red")

amp_genes <- out$amp;
sd_genes <- out$sigma;
phi_genes <- out$phi;

plot(density(phi_genes), col="red", main="Density plot of the phases")
plot(density(amp_genes), col="red", main="Density plot of the phases")
plot(density(sd_genes), col="red", main="Density plot of the phases")

ESS <- amp_genes^2; RSS <- sd_genes^2
SNR <- ESS/RSS;
plot(SNR, col="red", pch=20, lwd=1)

which(SNR > 500)
which.max(SNR)

plot(times, expr[, which.max(SNR)])
points(times, temp)
