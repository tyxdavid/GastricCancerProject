setwd("subcluster/Epithelial")

library(ggplot2)
library(Seurat)

source("Utils.R")


library(infercnv)

d<-readRDS("integrated_harmony.rds")

dd<-subset(d,Histological %in% c("normal","tumor","metastasis"))  ## only normal and tumor samples

rm(d)
gc()

print(table(dd[["Histological"]]))

## CNA score calculated as the original article
## InferCNV prep
dir.create("InferCNV")
setwd("InferCNV")

exp_matrix<-GetAssayData(object = dd, assay="RNA", slot = "counts")

annot<-dd@meta.data
annot$cell<-Cells(dd)

annot$Col2<-paste(annot$main_cluster,annot$Histological,annot$Patient,sep="_")
write.table(annot[,c("cell","Col2")],"annot_InferCNV_NormalRef.tsv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
ref_group<-unique(annot$Col2[which((annot$Histological %in% c("normal")))])  ## use cells from normal samples as reference



CNVObj<-CreateInfercnvObject(exp_matrix,"gen_pos.txt","annot_InferCNV_NormalRef.tsv",ref_group)
CNVObj <- run(CNVObj,
              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
              out_dir="InferCNV_NormalRef",  # dir is auto-created for storing outputs
              cluster_by_groups=FALSE,
              denoise=TRUE,
              noise_filter=0.05,
              HMM=FALSE,
              HMM_type="i3",
              remove_genes_at_chr_ends=TRUE,
              no_plot=TRUE,
              num_threads=20)
gc()


saveRDS(CNVObj,file="InferCNV_NormalRef_Epithelial.rds")
