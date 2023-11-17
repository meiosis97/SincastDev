testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- AddSincastObject(testdata)
testobj@misc
GetSincastObject(testobj)
GetSincastAssays(testobj)
CheckSincastObject(testobj)
CleanSincastAssays(testobj)


GetSincastAssays(testobj, assay = 'pseudobulk') <- testdata
CleanSincastAssays(testdata)@misc
CleanSincastAssays(testobj)@misc
GetSincastAssays(testobj, assay = 'pseudobulk')
testobj <- SincastAggregate(testobj, replace = T)
GetSincastAssays(testobj, assay = 'pseudobulk')



require(Sincast)
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testdata <- as.SincastSeurat(testdata)
testdata
