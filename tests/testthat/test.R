testdata <- readRDS(file = 'C:/Users/yidid/SincastDev/data/testdata.rds')
testobj <- AddSincast(testdata)
testobj@misc$SincastAssays@pseudobulk <- list(a = 1)
testobj <- AddSincast(testobj)

