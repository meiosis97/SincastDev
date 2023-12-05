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
testobj <- SincastImpute(testobj)
testobj@SincastAssays@pseudobulk@misc
testobj[["imputation"]]
testobj[["pseudobulk"]] <- NULL
testobj
