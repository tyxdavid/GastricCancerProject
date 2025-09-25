setwd('./main_cluster/subcluster/Epithelial')

library(Seurat)
library(CytoTRACE)
library(ggplot2)

d<-readRDS('seurat_for_monocle.rds')

exps<-GetAssayData(d,slot='count')
exps<-exps[rowSums(exps)>0,]

exps<-as.matrix(exps)

results <- CytoTRACE(exps, ncores = 48, enableFast=TRUE,subsamplesize=3000)

saveRDS(results,'CytoTRACE_results.rds')

res<-data.frame(CCAT=results$CytoTRACE,subcluster=d[['subcluster']]$subcluster)

uq_subcluster<-as.character(unique(d[['subcluster']]$subcluster))
medians<-numeric(length(uq_subcluster))
names(medians)<-uq_subcluster

for (i in uq_subcluster){
  medians[i]<-median(res$CCAT[res$subcluster==i])
}

medians<-sort(medians,decreasing=TRUE)
res$subcluster<-factor(res$subcluster,names(medians),ordered=TRUE)

g<-ggplot(res)+
   geom_boxplot(aes(x=subcluster,y=CCAT,fill=subcluster))+
   theme(legend.position='none',axis.text.x=element_text(angle=90,vjust=0.3))+
   scale_fill_hue()+ylab('CytoTRACE score')

ggsave(plot=g,'CytoTRACE_clusters.png')


