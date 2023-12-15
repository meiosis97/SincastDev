library("Sincast")
require(Matrix)
require(dplyr)
require(ggplot2)
require(Seurat)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj, replace = T)

p1 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = testdata$RNA_snn_res.0.15625)) + scale_color_discrete("")
p2 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = testdata$RNA_snn_res.0.883883476483184)) + scale_color_discrete("")
p3 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = testdata$RNA_snn_res.2.43569644814373)) + scale_color_discrete("")

p1 + p2 +p3




library(RANN)
library(glmnet)
nn.idx <- nn2(testdata@reductions$pca@cell.embeddings, k = 50)$nn.idx
b <- list()
response.gene <- "EGR2"
predictor.genes <- c("VCAN", "EGR1", "CCL3", "EGR3","CLEC9A",
                     "CCL4", "S100A12", "TNF", "CCL2",
                     "CCL5", "CD34", "CD36")

Y <- testobj[["imputation"]]@assays$RNA$counts[response.gene,] %>% scale()
X <- t(testobj[["imputation"]]@assays$RNA$counts[predictor.genes,])  %>% scale()
for(i in 1:nrow(nn.idx)){
  idx <- nn.idx[i,]
  y <- Y[idx]
  x <- X[idx,]
  b[[i]] <- lm(y~x)$coefficient
}
b <- t(do.call(rbind,b))
b <- b[-1,]
rownames(b) <- predictor.genes
colnames(b) <- colnames(testdata)

p1 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b["CCL3",])) +
  scale_color_gradientn("CCL3", colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p2 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b["EGR3",])) +
  scale_color_gradientn("EGR3", colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
p1+p2



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

b.imputed <- impute.function(b, 3)
limits <- quantile(b.imputed, c(0,1))
limits <- c(-1,1) * max(abs(limits))

b.imputed.residual <- t(apply(b.imputed, 1, function(x) (x-mean(x))/abs(mean(x))))
b.imputed.residual <- sign(b.imputed.residual) * abs(b.imputed.residual) ^ 0.25


cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"Spectral")))(100)

gene <- "CD83"
b.gene <-  b.imputed[gene,]
b.gene.residual <- b.imputed.residual[gene,]
limits <- c(-1,1) * max(abs(b.gene))
limits.residual <- c(-1,1) * max(abs(b.gene.residual))

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b.gene)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits)

ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = b.gene.residual)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")),
                        limits = limits.residual)


mu <- apply(b.imputed, 1, function(x) mean(abs(x)))
cv <- apply(b.imputed, 1, function(x) sd(x)/mean(x))
mu.cv.df <- data.frame(gene = names(mu),mu = mu, log.abs.cv = log(abs(cv)))

ggplot() + geom_point(data = mu.cv.df, aes(mu, log.abs.cv)) +
  geom_point(data = subset(mu.cv.df, log.abs.cv<0), aes(mu, log.abs.cv), col = "red")  +
  geom_text(data = subset(mu.cv.df, log.abs.cv<0), aes(mu, log.abs.cv, label = gene) ,col = "red")+
  scale_y_reverse() +
  geom_hline(yintercept = 0, linetype = "dashed") + scale_x_log10()


on <- apply(b,2,function(x) rank(-x) < 20)
on.imputed <- impute.function(on, 3)

p1 <- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = on["EGR3",])) +
  scale_color_manual(values = c("lightgrey","black"))
p2<- ggplot(data = as.data.frame(testdata@reductions$umap@cell.embeddings)) +
  geom_point(aes(umap_1, umap_2, col = on.imputed["EGR3",])) +
  scale_color_gradient(low = "lightgrey", high = "black")
p1+p2


pca.mod <- prcomp(t(on.imputed))
ggplot(data = as.data.frame(pca.mod$x)) +
  geom_point(aes(PC1, PC2, col = testdata$seurat_clusters))+scale_color_discrete("")

testdata
testobj <- as.Sincast(testdata)
testobj <- SincastImpute(testobj)
testobj <- SincastAggregate(testobj, size.factor = 10)
testobj@SincastAssays@original

testobj[["original"]]
ImputationPlot(testobj)
