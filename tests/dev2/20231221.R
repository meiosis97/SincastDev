library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
require(Seurat)
require(RANN)
require(glmnet)
require(randomForest)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
rownames(testdata)[which(testdata@assays$RNA@meta.data$var.features.rank < 100)]


my.pca <- function(x, y = NULL, k = 100){
  x <- scale(x)
  svd.x <- RSpectra::svds(x, k = k)
  s <- abs(diff(svd.x$d))
  mu <- mean(s[(0.8*k):k-1])
  sigma <- sd(s[(0.8*k):k-1])
  sk <- mu + 6*sigma
  rank <- which(s < sk)[1]
  if(length(rank)==0){
    warning(rank, " single values could be insufficient to approximate the data.")
    rank <- k
  }

  pcs <- x %*% svd.x$v[,1:rank]
  loadings <- svd.x$v[,1:rank]
  colnames(loadings) <- colnames(pcs) <- paste("PC", 1:rank, sep ='')
  rownames(pcs) <- rownames(x)
  rownames(loadings) <- colnames(x)

  lra <- tcrossprod(pcs, loadings)[,y]

  list(lra = lra, pcs = pcs, loadings = loadings)
}
pca.resust <- my.pca(t(testdata@assays$RNA$scale.data), y = "ISG15")


Y <- pca.resust$lra
X <- pca.resust$pcs
b <- list()
nn.idx <- nn2(pca.resust$pcs, k = 100)$nn.idx
for(i in 1:nrow(nn.idx)){
  idx <- sample(nn.idx[i,],50, replace = T)
  y <- Y[idx]
  x <- X[idx,]
  b[[i]] <- glmnet(x, y, standardize = FALSE, standardize.response = FALSE)$beta[,20]
  print(i)
}
b <- t(do.call(rbind,b))
b.all <- pca.resust$loadings %*% b
rownames(b.all) <- rownames(testdata@assays$RNA$scale.data)

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b.all["MX1",])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))


sort(-rowMeans(abs(b.all)))[1:50]
ggplot(data.frame(t(testdata@assays$RNA$scale.data))) + geom_point(aes(CD74, Y))

plot(Y, )

