library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
require(Seurat)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)

impute.function <- function(x, t){
  for(i in 1:t) x <- x %*% testobj[["imputation"]]@graphs$RNA_Sincast
  x
}
response.gene <- "EGR2"
predictor.genes <- VariableFeatures(testdata)
predictor.genes <- predictor.genes[predictor.genes != "EGR2"]
nn.idx <- nn2(testdata@reductions$pca@cell.embeddings, k = 50)$nn.idx
b <- list()
Y <- testobj[["imputation"]]@assays$RNA$counts[response.gene,] %>% scale()
X <- t(testobj[["imputation"]]@assays$RNA$counts[predictor.genes,])  %>% scale()
for(i in 1:nrow(nn.idx)){
  idx <- sample(nn.idx[i,],25, replace = T)
  y <- Y[idx]
  x <- X[idx,]
  b[[i]] <- glmnet(x, y, standardize = FALSE)$beta[,50]
  print(i)
}
b <- t(do.call(rbind,b))
