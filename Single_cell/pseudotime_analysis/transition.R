setwd('./main_cluster/subcluster/Epithelial')

library(Seurat)
library(monocle)
source('./Utils.R')

HSMM<-readRDS("tmp.rds")

gc()

clustering_DEG_genes<-readRDS("all_DEG_genes_hist.rds")


  HSMM_ordering_genes <-row.names(clustering_DEG_genes)[1:5000]

  HSMM <-  setOrderingFilter(HSMM,ordering_genes = HSMM_ordering_genes)
  gc()

  HSMM <-reduceDimension(HSMM, method = 'DDRTree')
  gc()

  saveRDS(HSMM,"ordered_HSMM.rds")

  HSMM <- orderCells(HSMM)

  saveRDS(HSMM,"ordered_HSMM.rds")

rm(HSMM)
gc()

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    Root_counts <- table(pData(cds)$State, pData(cds)$Histological)[,"tumor"]
	Root_counts_met <- table(pData(cds)$State, pData(cds)$Histological)[,"metastasis"]
	Root_counts_all <- table(pData(cds)$State)
	Root_counts_prop<-Root_counts/Root_counts_met
	names(Root_counts_prop)<-names(Root_counts)
	print(Root_counts_prop)
	if (anyNA(Root_counts_prop)){
	  Root_counts_prop<-Root_counts_prop[-which(is.na(Root_counts_prop))]
	}
	print(Root_counts_prop)
    return(as.numeric(names(Root_counts_prop)[which
          (Root_counts_prop == max(Root_counts_prop))]))
  } else {
    return (1)
  }
}


HSMM<-readRDS("ordered_HSMM.rds")

HSMM <- orderCells(HSMM, root_state = GM_state(HSMM))
gc()

saveRDS(HSMM,"ordered_HSMM.rds")
