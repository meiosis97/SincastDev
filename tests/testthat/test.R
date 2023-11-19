library("Sincast")
testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- as.SincastSeurat(testdata)
testobj@misc
GetSincastObject(testobj)
GetSincastAssays(testobj)
CheckSincastObject(testobj)
CleanSincastAssays(testobj)


GetSincastAssays(testobj, assay = 'pseudobulk') <- testdata
CleanSincastAssays(testdata)@Sincast
CleanSincastAssays(testobj)@Sincast
CheckSincastObject(testobj, assay = 'pseudobulk')
GetSincastObject(testdata)
GetSincastAssays(testobj)
testobj <- SincastAggregate(testobj, replace = T)
GetSincastAssays(testobj, assay = 'pseudobulk')@misc




testobj

