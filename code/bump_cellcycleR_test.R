
############  bump cell cycle example  ##################

data <- get(load("../output/cellcyler_data_liver_mash.rda"))

times <- data$times
expr <- data$data

celltime_levels=100
num_iter = 100

cycle_data <- expr
G <- dim(cycle_data)[2];
numcells <- dim(cycle_data)[1];
sigma <- array(0,G);
amp <- array(0,G);
phi <- array(0,G);

cell_times_iter <- times

# Fit linear models for each gene $g$ given the cell times [ linear model depends on fix.phase]

lmfit_list <- lapply(1:G, function(g)
{
  temp1 <- scales::rescale(cycle_data[,g], to=c(0,1))
  fit1 <- lm(temp1  ~ sin(0.5*cell_times_iter) + cos(0.5*cell_times_iter) -1);
  plot(cell_times_iter, fit1$fitted.values)
  points(cell_times_iter, temp1, col="blue")

  temp2 <- scales::rescale(cycle_data[,g], to = c(-1, 0))
  fit2 <- lm(temp2  ~ sin(0.5*cell_times_iter) + cos(0.5*cell_times_iter) -1);
  plot(cell_times_iter, fit2$fitted.values)
  points(cell_times_iter, temp2, col="blue")
  s1 <- summary(fit1)$r.squared
  s2 <- summary(fit2)$r.squared
  if( s1 > s2 ){
    fit <- fit1
    scale <- 0
  }else{
    fit <- fit2
    scale <- 1
  }

  out_sigma <- sd(fit$residuals);
  beta1 <- fit$coefficients[1];
  beta2 <- fit$coefficients[2];
  if(beta1==0 & beta2==0){
    stop(paste0("You have a gene with all 0 counts at gene",g));
  }
  out_amp <- sqrt(beta1^2 + beta2^2);
  out_phi <- atan3(as.numeric(beta2), as.numeric(beta1));
  ll <- list("out_amp"=out_amp, "out_phi"=out_phi, "out_sigma"=out_sigma,
             "out_scale" = scale)
  return(ll)
})

amp <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_amp))));
phi <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_phi))));
sigma <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_sigma))));
scale <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_scale))));

plot(times, expr[,467])

cycle_data_scaled <- matrix(0, dim(cycle_data)[1], dim(cycle_data)[2])
for(l in 1:dim(cycle_data_scaled)[2]){
  if(scale[l] == 0){
    cycle_data_scaled[,l] <- scales::rescale(cycle_data[,l], to=c(0,1))
  }else{
    cycle_data_scaled[,l] <- scales::rescale(cycle_data[,l], to=c(-1,0))
  }
}

plot(times, amp[300]*sin(0.5*cell_times_iter + phi[300]))
points(times, cycle_data_scaled[,300], col="blue", pch=20)

cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
num_celltime_class <- length(cell_times_class);

sin_class_times <- sin(0.5*cell_times_class);
cos_class_times <- cos(0.5*cell_times_class);
sin_phi_genes <- sin(phi);
cos_phi_genes <- cos(phi);
sinu_signal <- cbind(sin_class_times, cos_class_times) %*% rbind(amp*cos_phi_genes, amp*sin_phi_genes);
options(digits=12)
signal_intensity_per_class <- matrix(0, numcells, num_celltime_class)

signal_intensity_per_class <- do.call(rbind,lapply(1:numcells, function(cell)
{
  res_error <- sweep(sinu_signal,2,cycle_data_scaled[cell,]);
  res_error_adjusted <- -(res_error^2);
  res_error_adjusted <- sweep(res_error_adjusted, 2, 2*sigma^2, '/');
  out <- rowSums(sweep(res_error_adjusted,2,log(sigma)) - 0.5*log(2*pi));
  return(out)
}));

signal_intensity_class_exp <- do.call(rbind,lapply(1:dim(signal_intensity_per_class)[1], function(x)
{
  out <- exp(signal_intensity_per_class[x,]- max(signal_intensity_per_class[x,]));
  return(out)
}));

plot(cell_times_class, signal_intensity_per_class[5,])
abline(v = cell_times_iter[5], col="blue")

plot(cell_times_class, signal_intensity_class_exp[5,])
abline(v = cell_times_iter[5], col="blue")



data <- get(load("../output/cellcyler_data_liver_mash.rda"))

times <- data$times
expr <- data$data

out <- bump_cell_ordering_class(expr, celltime_levels = 100,
                                num_iter=100,
                                start = NULL,
                                verbose = TRUE,
                                save_path="../output/cell_order_liver_bump.rda")


plot(times, out$cell_times)
plot(out$cell_times, expr[,10])
plot(times, expr[,10], col="red")

cycle_data <- expr

celltime_levels=100
num_iter = 100

iter <- 1

