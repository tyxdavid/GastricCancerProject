setwd("./subcluster/Epithelial/InferCNV")
dir.create("CNV_analysis")
setwd("CNV_analysis")

library(Seurat)
library(ggplot2)
source("Utils.R")

d<-readRDS("./subcluster/Epithelial/integrated_harmony.rds")

dd<-subset(d,Histological %in% c("normal","tumor","metastasis"))

cnv_res<-readRDS("./subcluster/Epithelial/InferCNV/InferCNV_NormalRef_Epithelial.rds")
cnv_res<-cnv_res@expr.data
cnv_res<-cnv_res[,which(colnames(cnv_res) %in% colnames(dd))]

denoise.range<-c(0.95,1.05)
cnv_res[cnv_res>denoise.range[1] & cnv_res<denoise.range[2]]<-1

cnv_assay <- CreateAssayObject(data = cnv_res)

dd[["CNV"]]<-cnv_assay

DefaultAssay(dd) <- "CNV"


saveRDS(dd,"CNV_analysis.rds")
