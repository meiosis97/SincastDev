testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testdata


testobj <- CreateSincastAssays(testdata)
Seurat::Misc(testobj, slot = "SincastAssays") <- list(a = 1,b=2,c=3)
testobj <- CreateSincastAssays(testobj, replace = T)
testobj <- SincastAggregate(testobj, size.factor = 10, replace = T)


testobj@misc
testobj@misc$SincastAssays@pseudobulk[[1]]
validObject(testobj@misc$SincastAssays, test = T)


library(Sincast)
names(testobj@misc)
