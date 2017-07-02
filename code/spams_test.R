
###########  spams testing  #################

library(spams)

library(spams)
library(png)
I = readPNG(paste("../../spams/extdata/boat.png"))

m = 8;n = 8;
X = spams.im2col_sliding(I,m,n)

X = X - matrix(rep(colMeans(X),nrow(X)),nrow(X),ncol(X),byrow = T)
X = X / matrix(rep(sqrt(colSums(X*X)),nrow(X)),nrow(X),ncol(X),byrow=T)
lambda1 = 0.15

tic = proc.time()
D <- spams.trainDL(X,K = 100,lambda1 = lambda1, numThreads = 1, batchsize = 400,
                   iter = 1000)
tac = proc.time()
t = (tac - tic)[['elapsed']]
.printf("time of computation for Dictionary Learning: %f\n",t)
.objective(X,D,lambda1)

