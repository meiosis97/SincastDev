library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
require(Seurat)
require(RANN)
require(glmnet)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)


impute.function <- function(x, t){
  for(i in 1:t) x <- x %*% t(testobj["imputation"]@graphs$sincast_RNA_DO)
  x
}
response.gene <- "EGR2"
predictor.genes <- VariableFeatures(testdata)
predictor.genes <- predictor.genes[predictor.genes != "EGR2"]
nn.idx <- nn2(testdata@reductions$pca@cell.embeddings, k = 100)$nn.idx
b <- list()
Y <- testobj["imputation"]@assays$RNA$counts[response.gene,] %>% scale()
X <- t(testobj["imputation"]@assays$RNA$counts[predictor.genes,])  %>% scale()
for(i in 1:nrow(nn.idx)){
  idx <- sample(nn.idx[i,],25, replace = T)
  y <- Y[idx]
  x <- X[idx,]
  b[[i]] <- glmnet(x, y, standardize = FALSE, standardize.response = FALSE)$beta[,50]
  print(i)
}
b <- t(do.call(rbind,b))
b.imputed <- impute.function(b, 3)
b.imputed <- LowRankApprox(abs(b.imputed))
b.imputed <- sign(b) * b.imputed


gene <- "CCL4"
b.gene <-  b.imputed[gene,]
limits <- c(-1,1) * max(abs(b.gene))
limits <- c(-1,1) * max(abs(b.gene))

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b.gene)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits)



