library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)




x <- t(testdata@assays$RNA$data)
test <- RSpectra::svds(x, 100)
Xhat <- tcrossprod(test$u %*% diag(test$d), test$v)
tot.var <- sum(test$d)/sum(x^2)*nrow(x)
colnames(Xhat) <- colnames(x)

lambda1 <- (Xhat - x)^2
lambda2 <- (Xhat - t(testobj[["imputation"]]@assays$RNA$counts))^2
lambda3 <- lambda1/(lambda1+lambda2)
imp.data <-  x * (1-lambda3) + t(testobj[["imputation"]]@assays$RNA$counts) * lambda3

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = imp.data[,"EGR2"])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))


plot(t(testobj[["imputation"]]@assays$RNA$counts)[,"CCL3"], x[,"CCL3"])
plot(imp.data[,"CCL3"], x[,"CCL3"])
plot(Xhat[,"CCL3"], x[,"CCL3"])


plot(imp.data[,"CCL3"], imp.data[,"CCL4"])
plot(lambda3[,10])

Xres <- x - Xhat
test <- RSpectra::svds(Xres, 1)
Xhat <- Xhat + tcrossprod(test$u, test$v) * test$d


plot(Xhat[,"CCL3"],x[,"CCL3"])
