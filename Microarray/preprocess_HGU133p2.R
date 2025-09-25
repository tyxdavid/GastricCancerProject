library(affy)
library(hgu133plus2.db)
library(hgu133plus2cdf)

files<-list.files("./data/",pattern="gz")
sampleAccession<-substr(files,1,9)

d<-ReadAffy(filenames = paste0("./data/",files),compress = T,sampleNames = sampleAccession)

eset <- expresso(d, normalize.method="quantiles",
                 bgcorrect.method="rma",pmcorrect.method="pmonly",
                 summary.method="medianpolish",verbose = T)

d<-exprs(eset)

saveRDS(eset,'eset.rds')
saveRDS(d,'exps.rds')
