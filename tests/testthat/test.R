testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- AddSincastObject(testdata)
testobj@misc
GetSincastObject(testobj)
CheckSincastObject(testobj)
CleanSincastAssays(testobj)


GetSincastAssays(testobj, assay = 'pseudobulk') <- testdata
CleanSincastAssays(testdata)@misc
CleanSincastAssays(testobj)@misc
