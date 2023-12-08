library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)

library(RANN)
library(glmnet)
impute.function <- function(x, t){
  for(i in 1:t) x <- x %*% testobj[["imputation"]]@graphs$RNA_Sincast
  x
}

nn.idx <- nn2(testdata@reductions$pca@cell.embeddings, k = 50)$nn.idx
b <- list()
Y <- testobj[["imputation"]]@assays$RNA$counts["EGR2",] %>% scale()
X <- t(testobj[["imputation"]]@assays$RNA$counts[c("VCAN", "EGR1", "CCL3", "EGR3","CLEC9A",
                                                   "CCL4", "S100A12", "TNF", "CCL2",
                                                   "CCL5", "CD34", "CD36"),])  %>% scale()
for(i in 1:nrow(nn.idx)){
  idx <- sample(nn.idx[i,],25, replace = T)
  y <- Y[idx]
  x <- X[idx,]
  b[[i]] <- glmnet(x, y, alpha = 0, standardize = FALSE)$beta[,100]
}
b <- t(do.call(rbind,b))

b.imputed <- impute.function(b, 3)
limits <- quantile(b.imputed, c(0,1))
limits <- c(-1,1) * max(abs(limits))

b.imputed.residual <- t(apply(b.imputed, 1, function(x) (x-mean(x))/abs(mean(x))))
b.imputed.residual <- sign(b.imputed.residual) * abs(b.imputed.residual) ^ 0.25
limits.residual <- quantile(b.imputed.residual, c(0,1))
limits.residual  <- c(-1,1) * max(abs(limits.residual))
breaks.residual <- seq(limits.residual[1], limits.residual[2], length.out = 5)
labels.residual <- round(sign(breaks.residual) * breaks.residual^4,2)

cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"Spectral")))(100)

gene <- "CCL3"
b.gene <-  b.imputed[gene,]
b.gene.residual <- b.imputed.residual[gene,]

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b.gene)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits)

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b.gene.residual)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits.residual,
                        breaks = breaks.residual,
                        labels = labels.residual)


ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = testdata@assays$RNA$data[gene,])) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))

ggplot() +
  geom_point(aes(X[,gene], Y, col = testdata$seurat_clusters))



mu <- apply(b.imputed, 1, function(x) mean(abs(x)))
cv <- apply(b.imputed, 1, function(x) sd(x)/mean(x))
ggplot() +
  geom_text(aes(mu, abs(cv), label = names(mu))) + scale_y_reverse() +
  geom_hline(yintercept = 1, linetype = "dashed")



pca.mod <- prcomp(t(b.imputed))
ggplot(data = as.data.frame(pca.mod$x)) +
  geom_point(aes(PC1, PC2, col = testdata$seurat_clusters))
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits)

