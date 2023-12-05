library("Sincast")
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
CheckSincastObject(testobj)
CheckSincastAssays(testobj@SincastAssays)
GetSincastAssays(testobj, "original")
GetSincastAssays(testobj)
CleanSincastAssays(testobj@SincastAssays)
CleanSincastAssays(testobj)
GetSincastAssays(testobj, assay = 'pseudobulk') <- testdata
testobj <- SincastAggregate(testobj, replace = T, size.factor = 10)
testobj <- SincastImpute(testobj, replace = T)
testobj@SincastAssays@pseudobulk@misc
ImputationPlot(testobj, color.by = "EGR2", anno.by = c("ident", "CD34"), dims = 1:3)


################################################################################
library(RANN)
nn.idx <- nn2(testdata@reductions$pca@cell.embeddings, k = 60)$nn.idx
r <- c()
Y <- testobj[["imputation"]]@assays$RNA$counts["EGR2",]
X <- testobj[["imputation"]]@assays$RNA$counts["ERBB2",]
for(i in 1:nrow(nn.idx)){
  idx <- sample(nn.idx[i,],30)
  y <- scale(Y[idx])
  x <- scale(X[idx])
  r[i] <- lm(y~x)$coefficient[2]
}

r <- r
r.imputed <- impute.function(r, 3)
limits <- range(r.imputed)
limits <- c(-1,1) * max(abs(limits))
ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = r.imputed)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits)

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = X)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(X, Y, col = testdata$seurat_clusters))

impute.function <- function(x, t){
  out <- x
  for(i in 1:t) out <- out %*% testobj[["imputation"]]@graphs$RNA_Sincast
  mod <- lm(sort(x)~ sort(out))
  out <- as.numeric(out * coef(mod)[2] + coef(mod)[1])
  out
}

