library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)

# Test average
x <- t(testdata@assays$RNA$data)
svd.out <- RSpectra::svds(x, 20)
x.lra <- tcrossprod(svd.out$u %*% diag(svd.out$d), svd.out$v)
min.lra <- apply(x.lra, 2, min)
min.lra[min.lra>0] <- 0
x.lra <- sweep(x.lra, 2, min.lra)
colnames(x.lra) <- colnames(x)
QQScale <- function(before, after) {
  #QQ regression
  for(i in 1:nrow(before)){
    y <- before[i,]
    x <- after[i,y>0]
    y <- before[i,y>0]
    if(length(y)==1){
      after[i,] <- after[i,]*y/x
    }else if(length(y)>1){
      mod <- lm(sort(y)~0+sort(x))
      after[i,] <- as.numeric(after[i,] * coef(mod)[1])
    }
  }

  after
}
x.lra <- t(QQScale(t(x) %>% as.matrix(), t(x.lra)))


v1 <- (x-t(testobj[["imputation"]]@assays$RNA$counts))^2
v2 <- (x-x.lra)^2
x.bar <- (v2/(v1+v2)) * t(testobj[["imputation"]]@assays$RNA$counts) + (v1/(v1+v2)) * x.lra
v1 <- (x-x.bar)^2    #((x - x.lra)^2 + (x-t(testobj[["imputation"]]@assays$RNA$counts))^2)/2
v2 <- ((x.bar - x.lra)^2 + (x.bar-t(testobj[["imputation"]]@assays$RNA$counts))^2)/2
x.bar <- (v2/(v1+v2)) * x + (v1/(v1+v2)) * x.bar

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.bar[,3])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))


plot(x.bar[,"CCL3"], x.bar[,"CCL4"])

abline(0,1)
plot(x[,"CLEC9A"], t(testobj[["imputation"]]@assays$RNA$counts)[,"CLEC9A"])
abline(0,1)
plot(x[,3],  t(testobj[["imputation"]]@assays$RNA$counts)[,3])
abline(0,1)
plot(x[,3], x.bar[,3])
abline(0,1)

# Test low rank approximation
x <- t(testdata@assays$RNA$data)
svd.out <- RSpectra::svds(t(testobj[["imputation"]]@assays$RNA$counts), 30)
x.lra <- tcrossprod(svd.out$u %*% diag(svd.out$d), svd.out$v)
colnames(x.lra) <- rownames(testdata@assays$RNA$data)
min.lra <- apply(x.lra, 2, min)
x.lra <- sweep(x.lra, 2, min.lra)
colnames(x.lra) <- colnames(x)
QQScale <- function(before, after) {
  #QQ regression
  for(i in 1:nrow(before)){
    y <- before[i,]
    x <- after[i,y>0]
    y <- before[i,y>0]
    if(length(y)==1){
      after[i,] <- after[i,]*y/x
    }else if(length(y)>1){
      mod <- lm(sort(y)~0+sort(x))
      after[i,] <- as.numeric(after[i,] * coef(mod)[1])
    }
  }

  after
}
x.lra <- t(QQScale(t(x) %>% as.matrix(), t(x.lra)))


ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.lra[,"CCL3"])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))

plot(x.lra[,"CCL3"], x.lra[,"CCL4"])
abline(0,1)

ggplot() +
  geom_point(aes(x.lra[,"CCL3"], x.lra[,"CCL4"], col = testdata$seurat_clusters))
