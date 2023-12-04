library("Sincast")
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
CheckSincastObject(testobj)
CheckSincastAssays(testobj@SincastAssays)
GetSincastAssays(testobj, "original")
GetSincastAssays(testobj)
CleanSincastAssays(testobj)
GetSincastAssays(testobj, assay = 'pseudobulk') <- testdata
testobj <- SincastAggregate(testobj, replace = T, size.factor = 10)
testobj <- SincastImpute(testobj, replace = T)
testobj@SincastAssays@pseudobulk@misc
ImputationPlot(testobj, color.by = "EGR2", anno.by = c("ident", "CD34"), dims = 1:3)

Seurat::DimPlot(testdata)
Seurat::FeaturePlot(testdata, features = "CCL3")




################################################################################
library(RANN)
nn.idx <- nn2(testdata@reductions$pca@cell.embeddings, k = 30)$nn.idx
r <- c()
Y <- testobj[["imputation"]]@assays$RNA$counts["CCL3",]
X <- testobj[["imputation"]]@assays$RNA$counts["CCL4",]
for(i in 1:nrow(nn.idx)){
  idx <- nn.idx[i,]
  y <- Y[idx]
  x <- X[idx]
  r[i] <- cor(y,x)
}

r <- r
r.imputed <- impute.function(r, 3)

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = r.imputed)) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(11,"Spectral"), limits = c(-1, 1))


impute.function <- function(x, t){
  out <- x
  for(i in 1:t) out <- out %*% testobj[["imputation"]]@graphs$RNA_Sincast
  mod <- lm(sort(x)~ sort(out))
  out <- as.numeric(out * coef(mod)[2] + coef(mod)[1])
  out[out>1] <- 1
  out[out<-1] <- -1
  out
}

plot(X,Y)
