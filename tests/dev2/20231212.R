library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
require(plotly)
require(ggpubr)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T, do.post.scale = F, preserve.zero = F)
testobj <- SincastAggregate(testobj, replace = T, size.factor = 10,n.pool = 30)
x <- t(testdata@assays$RNA$data)
x.imp <- t(testobj[["imputation"]]@assays$RNA$counts)
rankTrans <- function(data){
  (apply(data,2,rank, ties.method = 'min')-1)/(nrow(data) - 1)
}
pca <- prcomp(t(rankTrans(testobj[["pseudobulk"]]@assays$RNA$counts)))


# Original data
predicted <- predict(pca, t(rankTrans(t(x))))
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,
            marker = list(symbol = "x", size = 5))
p1 <- ggplot() + geom_point(aes(x[,"CCL3"], x[,"CCL4"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("CCL3") + ylab("CCL4") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))
p2 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x[,1])) +
  scale_color_gradientn("AL627309.1",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p3 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x[,"EGR2"])) +
  scale_color_gradientn("EGR2",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
ggarrange(NULL,p1,p2,p3, align = c('hv'))


# Low rank approximation
LowRankApprox <- function(x){
  svd.x <- RSpectra::svds(x, 100)
  s <- abs(diff(svd.x$d))
  mu <- mean(s[79:99])
  sigma <- sd(s[79:99])
  sk <- mu + 6*sigma
  k <- which.min(s > sk)
  x.lra <- tcrossprod(x %*% svd.x$v[,1:k], svd.x$v[,1:k])
  colnames(x.lra) <- colnames(x)
  x.lra
}
x.lra <- LowRankApprox(x)

predicted <- predict(pca, t(rankTrans(t(x.lra))))
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,
            marker = list(symbol = "x", size = 5))
p1 <- ggplot() + geom_point(aes(x.lra[,"CCL3"], x.lra[,"CCL4"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("CCL3") + ylab("CCL4") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))
p2 <- ggplot() + geom_point(aes(x.lra[,"CCL3"], x[,"CCL3"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("Imputed CCL3") + ylab("Original CCL3") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))+
  geom_abline(slope = 1, intercept = 0)
p3 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.lra[,1])) +
  scale_color_gradientn("AL627309.1",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p4 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.lra[,"EGR2"])) +
  scale_color_gradientn("EGR2",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
ggarrange(p1,p2,p3,p4, align = c('hv'))


# Sincast + Median scaling
qq.scale <- function(X, Y) {
  scale.factor <- apply(X, 1, function(x) median(x[x != 0])) /
      apply(replace(Y, X == 0, NA), 1, function(y) median(y, na.rm = T))

  Y <- Y * scale.factor
  Y[is.na(Y)] <- 0
  Y
}
x.imp.qq <-  t(qq.scale(t(x) %>% as.matrix, t(x.imp) %>% as.matrix))

predicted <- predict(pca, t(rankTrans(t(x.imp.qq))))
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,
            marker = list(symbol = "x", size = 5))
p1 <- ggplot() + geom_point(aes(x.imp.qq[,"CCL3"], x.imp.qq[,"CCL4"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("CCL3") + ylab("CCL4") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))
p2 <- ggplot() + geom_point(aes(x.imp.qq[,"CCL3"], x[,"CCL3"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("Imputed CCL3") + ylab("Original CCL3") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))+
  geom_abline(slope = 1, intercept = 0)
p3 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.imp.qq[,1])) +
  scale_color_gradientn("AL627309.1",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p4 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.imp.qq[,"EGR2"])) +
  scale_color_gradientn("EGR2",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
ggarrange(p1,p2,p3,p4, align = c('hv'))


# Sincast
predicted <- predict(pca, t(rankTrans(t(x.imp))))
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,
            marker = list(symbol = "x", size = 5))
p1 <- ggplot() + geom_point(aes(x.imp[,"CCL3"], x.imp[,"CCL4"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("CCL3") + ylab("CCL4") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))
p2 <- ggplot() + geom_point(aes(x.imp[,"CCL3"], x[,"CCL3"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("Imputed CCL3") + ylab("Original CCL3") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))+
  geom_abline(slope = 1, intercept = 0)
p3 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.imp[,1])) +
  scale_color_gradientn("AL627309.1",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p4 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.imp[,"EGR2"])) +
  scale_color_gradientn("EGR2",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
ggarrange(p1,p2,p3,p4, align = c('hv'))


# Zero preserving Low rank approximation
LowRankApprox <- function(x){
  svd.x <- RSpectra::svds(x, 100)
  s <- abs(diff(svd.x$d))
  mu <- mean(s[79:99])
  sigma <- sd(s[79:99])
  sk <- mu + 6*sigma
  k <- which.min(s > sk)
  x.lra <- tcrossprod(x %*% svd.x$v[,1:k], svd.x$v[,1:k])
  x.lra <- apply(x.lra, 2, function(x){
    x[x<abs(quantile(x[x<0], 0.1))] <- 0
    x
  } )
  colnames(x.lra) <- colnames(x)
  x.lra
}
x.lra.0 <- LowRankApprox(x)

predicted <- predict(pca, t(rankTrans(t(x.lra.0))))
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,
            marker = list(symbol = "x", size = 5))
p1 <- ggplot() + geom_point(aes(x.lra.0[,"CCL3"], x.lra.0[,"CCL4"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("CCL3") + ylab("CCL4") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))
p2 <- ggplot() + geom_point(aes(x.lra.0[,"CCL3"], x[,"CCL3"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("Imputed CCL3") + ylab("Original CCL3") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))+
  geom_abline(slope = 1, intercept = 0)
p3 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.lra.0[,1])) +
  scale_color_gradientn("AL627309.1",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p4 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.lra.0[,"EGR2"])) +
  scale_color_gradientn("EGR2",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
ggarrange(p1,p2,p3,p4, align = c('hv'))


# Zero preserving Low rank approximation on Sincast
LowRankApprox <- function(x){
  svd.x <- RSpectra::svds(x, 100)
  s <- abs(diff(svd.x$d))
  mu <- mean(s[79:99])
  sigma <- sd(s[79:99])
  sk <- mu + 6*sigma
  k <- which.min(s > sk)
  x.lra <- tcrossprod(x %*% svd.x$v[,1:k], svd.x$v[,1:k])
  x.lra <- apply(x.lra, 2, function(x){
    x[x<abs(quantile(x[x<0], 0.1))] <- 0
    x
  } )
  colnames(x.lra) <- colnames(x)
  x.lra
}
x.imp.0 <- LowRankApprox(x.imp)

predicted <- predict(pca, t(rankTrans(t(x.imp.0))))
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,
            marker = list(symbol = "x", size = 5))
p1 <- ggplot() + geom_point(aes(x.imp.0[,"CCL3"], x.imp.0[,"CCL4"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("CCL3") + ylab("CCL4") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))
p2 <- ggplot() + geom_point(aes(x.imp.0[,"CCL3"], x[,"CCL3"], col = testobj[["imputation"]]$seurat_clusters)) + xlab("Imputed CCL3") + ylab("Original CCL3") +
  scale_color_manual("cluster", values = RColorBrewer::brewer.pal(8, "Set2"))+
  geom_abline(slope = 1, intercept = 0)
p3 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.imp.0[,1])) +
  scale_color_gradientn("AL627309.1",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p4 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.imp.0[,"EGR2"])) +
  scale_color_gradientn("EGR2",colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
ggarrange(p1,p2,p3,p4, align = c('hv'))

