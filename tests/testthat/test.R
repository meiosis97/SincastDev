library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
require(Seurat)
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
testobj[["imputation"]]
testobj[["pseudobulk"]] <- NULL
testobj["pseudobulk"]
testobj
ImputationPlot(testobj, color.by = "EGR2", anno.by = c("ident", "CD34"), dims = 1:3)
testobj <- BuildSincastAtlas(testobj, replace = T)
summary(testobj)
testobj@SincastAtlas@pseudobulk@misc
capture.output(testobj@SincastAtlas)
testobj@summary@active.assay

AtlasPlot(testobj, atlas = "original")
help("AtlasPlot")

dim(testobj@SincastAtlas@original@assays$RNA$counts)
tes <- function(){
  message("\r  Now construct affinity matrix. ", appendLF = F)
  message(paste("\r  construct affinity matrix: Calcualting distance."), appendLF = F)
  message(paste("\r  Construct affinity matrix: Scaling distance."), appendLF = F)
}
tes()

testdata <- DietSeurat(testdata, layers = c('counts'),)
LayerData(testdata, layer = "data") <- NULL

testdata@assays$R

do.call(ScaleData, args)
ScaleData()
help(ScaleData)
rlang:::exec(ScaleData, !!!args)

rlang::exec(ScaleData, !!!args)

SeuratObject::LayerData(testobj[["imputation"]], layer = "data") <- SeuratObject::LayerData(testobj[["imputation"]], layer = "counts")

help("timestamp")
summary(testobj)
timestamp(quiet = F)
difftime(timestamp(quiet = F), timestamp(quiet = F))
################################################################################
library(RANN)
library(glmnet)
impute.function <- function(x, t){
  for(i in 1:t) x <- x %*% testobj[["imputation"]]@graphs$RNA_Sincast
  x
}

nn.idx <- nn2(testdata@reductions$pca@cell.embeddings, k = 50)$nn.idx
b <- list()
Y <- testobj[["imputation"]]@assays$RNA$counts["CCL3",] %>% scale()
X <- t(testobj[["imputation"]]@assays$RNA$counts[c("EGR2","CLEC9A","CCL4", "S100A12", "TNF", "CCL2", "CCL5", "CD34"),])  %>% scale()
for(i in 1:nrow(nn.idx)){
  idx <- sample(nn.idx[i,],25, replace = T)
  y <- Y[idx]
  x <- X[idx,]
  b[[i]] <- glmnet(x, y, alpha = 0, standardize = FALSE)$beta[,100]
}
b <- t(do.call(rbind,b))

b.imputed <- impute.function(b, 3)
limits <- range(b.imputed)
limits <- c(-1,1) * max(abs(limits))
b.imputed.residual <- t(apply(b.imputed, 1, function(x) (x-mean(x))/mean(x)))
limits.redisual <- range(b.imputed.residual)
limits.redisual  <- c(-1,1) * max(abs(limits.redisual))

gene <- "CCL5"
b.gene <-  b.imputed[gene,]
b.gene.residual <- b.imputed.residual[gene,]

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b.gene)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits)

# ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
#   geom_point(aes(umap_1, umap_2, col = b.gene.residual)) +
#   scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
#
# ggplot() +
#   geom_point(aes(gene, Y, col = testdata$seurat_clusters))

mu <- apply(b.imputed, 1, function(x) mean(abs(x)))
cv <- apply(b.imputed, 1, function(x) sd(x)/mean(x))
ggplot() +
  geom_text(aes(mu, 1/abs(cv), label = names(mu)))

################################################################################
setGeneric("test", function(a = 1, b = 2, c = 4, ...) {
  standardGeneric("test")
})

setMethod("test", "numeric", definition = function(a = 1, d= 5, b = 6, c = 8, ...){
  message(a,b,c,d)
})


test(a= 1)
