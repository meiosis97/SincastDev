library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
require(plotly)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)
testobj <- SincastAggregate(testobj, replace = T, size.factor = 10,n.pool = 30)

rankTrans <- function(data){
  (apply(data,2,rank, ties.method = 'min')-1)/(nrow(data) - 1)
}

pca<-prcomp(t(rankTrans(testobj[["pseudobulk"]]@assays$RNA$counts)))
predicted <- predict(pca, t(rankTrans(testobj[["imputation"]]@assays$RNA$counts)))
ggplot()+geom_point(data = data.frame(pca$x), aes(PC1,PC2, col = substr(testobj[["pseudobulk"]]$agg.label, 1,1)))+
  geom_point(data = data.frame(predicted), aes(PC1,PC2, col = testobj[["imputation"]]$seurat_clusters), shape = "x")
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,  marker = list(symbol = "x"))



testdata@reductions$pca@assay.used
Seurat::Reductions(testdata, "pca")@assay.used























x <- t(testdata@assays$RNA$data)
x.imp <- t(testobj[["imputation"]]@assays$RNA$counts)

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
  x.lra
}
qq.scale <- function(before, after){
  for(i in 1:nrow(before)){
    y <- before[i,]
    x <- after[i,]
    npos <- sum(y>0)
    if(npos>1){
      mod <-  lm(sort(y)~  0+sort(x))
      after[i,] <- as.numeric(after[i,] * coef(mod)[1])
    }
  }
  after
}

x.lra <- LowRankApprox(x)
x.imp.lra <- LowRankApprox(x.imp)
x.hat <- qq.scale(x.lra, x.imp.lra)
x.hat <- (x.hat + x.lra)/2
colnames(x.lra) <- colnames(x.hat) <- colnames(x.imp.lra) <- colnames(x)


svd.x.imp <- RSpectra::svds(x, 100)


RSpectra::



  mu <- apply(x,2,mean)
sigma <- apply(x,2,sd)
plot(log(mu),log(sigma))






x.lra <- tcrossprod(x.imp %*% svd.x$v, svd.x$v)
# x.postimp <- tcrossprod(x.imp %*% svd.x$v, svd.x$v)
# min.x.lar <- apply(x.lra, 2, min)
# min.x.lar[min.x.lar>0]<-0
# x.postimp <- sweep(x.postimp,2,min.x.lar)
# x.lra <- sweep(x.lra,2,min.x.lar)
# colnames(x.postimp) <- colnames(x)
colnames(x.lra) <- colnames(x)

# x.lra <- t(qq.scale(t(x) %>% as.matrix, t(x.lra) %>% as.matrix))
# x.imp <- t(qq.scale(t(x) %>% as.matrix, t(x.imp) %>% as.matrix))
x.lra <- apply(x.lra, 2, function(x){
  x[x<abs(quantile(x[x<0], 0.1))] <- 0
  x
} )

qq.scale <- function(before, after){
  for(i in 1:nrow(before)){
    y <- before[i,]
    x <- after[i,]
    npos <- sum(y>0)
    if(npos==1){
      after[i,] <- after[i,] * y/x
    }else if(npos>1){
      mod <-  lm(sort(y)~  0+sort(x))
      after[i,] <- as.numeric(after[i,] * coef(mod)[1])
    }
  }
  after
}

x.bar <- t(qq.scale(t(x.lra) %>% as.matrix, t(x.imp) %>% as.matrix))
x.bar <- x.bar * (x.lra!=0)



quantile(x[x<0], 0.1)



#
s <- abs(diff(svd.x$d))
mu <- mean(s[79:99])
sigma <- sd(s[79:99])
sk <- mu + 6*sigma
k <- which.min(s > sk)

#
ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.hat[,"CLEC9A"])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))


#
plot(x.imp[,"CCL3"], x.imp[,"CCL4"])
plot(x.imp[,"CCL3"], x.imp[,"EGR2"])

plot(x.lra[,"CCL3"], x.lra[,"CCL4"])
plot(x.bar[,"CCL3"], x.bar[,"EGR2"])


#
plot(x.imp.lra[,"CCL3"], x[,"CCL3"])
plot(x.imp[,"EGR2"], x[,"EGR2"])
plot(x.lra[,"CLEC9A"], x[,"CLEC9A"])
abline(0,1)

plot(x.imp[,"CCL3"], x[,"CCL3"])
plot(x.imp[,"EGR2"], x[,"EGR2"])
plot(x.hat[,"CLEC9A"], x[,"CLEC9A"])
abline(0,1)




rankTrans <- function(data){
  (apply(data,2,rank, ties.method = 'min')-1)/(nrow(data) - 1)
}



pca<-prcomp(t(rankTrans(testobj[["pseudobulk"]]@assays$RNA$counts)))
predicted <- predict(pca, t(rankTrans(t(x.lra))))
ggplot()+geom_point(data = data.frame(pca$x), aes(PC1,PC2, col = substr(testobj[["pseudobulk"]]$agg.label, 1,1)))+
  geom_point(data = data.frame(predicted), aes(PC1,PC2, col = testobj[["imputation"]]$seurat_clusters), shape = "x")
require(plotly)
plot_ly() %>% add_trace(data = data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color =  substr(testobj[["pseudobulk"]]$agg.label, 1,1)) %>%
  add_trace(data = data.frame(predicted), x = ~PC1, y = ~PC2, z = ~PC3, color = testobj[["imputation"]]$seurat_clusters,  marker = list(symbol = "x"))
