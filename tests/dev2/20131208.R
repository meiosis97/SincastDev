library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)

# Test PLS
x <- scale(t(testdata@assays$RNA$data))
x.center <- attr(x, "scaled:center")
x.scale <- attr(x, "scaled:scale")
x.imp <-  scale(t(testobj[["imputation"]]@assays$RNA$counts), center = x.center, scale = x.scale)
rm(testobj)

svd.out <- RSpectra::svds(x, 20)
pcs <- svd.out$u %>% scale()
pcs.imp <- x.imp %*% svd.out$v %>% scale()
x.lra <- tcrossprod(svd.out$u %*% diag(svd.out$d), svd.out$v)
x.lra <- scale(x.lra, center = -x.center, scale = 1/x.scale)


colnames(x.lra) <- colnames(x)

pls.svd.out <- RSpectra::svds(t(pcs) %*% pcs.imp, 20)
x.postscale <- pcs.imp %*% pls.svd.out$v %*% t(pls.svd.out$v) %*% t(svd.out$v)
colnames(x.postscale) <- colnames(x)
#x.postscale <- scale(x.postscale, center = -x.center, scale = 1/x.scale)

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.[,"CCL3"])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))


plot(pcs %*% pls.svd.out$u[,1], pcs.imp %*% pls.svd.out$v[,1])

# x.imp <- x.imp %*% svd.out$v %*% t(svd.out$v)
# colnames(x.imp) <- colnames(x)
# x.imp <- scale(x.imp, center = -x.center, scale = 1/x.scale)
# x <- scale(x, center = -x.center, scale = 1/x.scale)
# rm(testobj)

# ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
#   geom_point(aes(umap_1, umap_2, col = x[,"EGR2"])) +
#   scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
plot(x.lra[,1],x.imp[,1])






######################
x <- t(testdata@assays$RNA$data)
x.imp <-  t(testobj[["imputation"]]@assays$RNA$counts)
ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = x.imp[,"CCL4"])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))

ggplot() + geom_point(aes(x.imp[,"CCL3"], x[,"CCL3"]))

