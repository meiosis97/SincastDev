testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- CreateSincastAssays(testdata)
Seurat::Misc(testobj, slot = "SincastAssays") <- list(a = 1,b=2,c=3)
testobj <- CreateSincastAssays(testobj)
testobj <- SincastAggregate(testobj, size.factor = 100)



testobj@misc
testobj@misc$SincastAssays@pseudobulk[[1]]
validObject(testobj@misc$SincastAssays, test = T)


