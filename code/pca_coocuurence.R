

text <- c("I am speaking", "I was eating", "I was travelling")

mat <- rbind(c(0, 1, 2,  0, 0, 0), c(1, 0, 0,  1, 0, 0), 
             c(2, 0, 0,  0, 1, 1), 
             c(0, 1, 0,  0, 0, 0), c(0, 0, 1, 0, 0, 0),
             c(0, 0, 1, 0, 0, 0) )

rownames(mat) <- c("I", "am", "was", "speaking", "eating", "travelling")
colnames(mat) <- rownames(mat)

pca <- tsne::tsne(mat)
plot(pca[,1], pca[,2], col="red")
text(pca[,1], pca[,2], rownames(mat))



pca <- svd(mat)
x <- jitter(pca$v[,1], factor=2)
y <- jitter(pca$v[,2], factor=2)
plot(x, y, col="red", xlim = c(-1, 1), type="n", xlab = "SVD1", ylab = "SVD2")
text(x, y, rownames(mat), col = "red", cex = 1.5)
