##################
##### Fig S2 #####
##################

setwd("./subcluster/Epithelial")

library(Seurat)
library(ggplot2)
source("Utils.R")

d<-readRDS("integrated_harmony.rds")



dd<-subset(d,Histological %in% c("tumor","normal","metastasis"))

CT<-dd@meta.data
CT<-subset(CT,Histological %in% c("normal"))


deg_malignant_res<-DE_malignant_method(dd,normal_cells=rownames(CT),num_markers=NA,log2FC=1,tmpdir="DE_malignant_log2FC1",
									   cell_tolerate=Inf)

saveRDS(deg_malignant_res,"DE_malignant_log2FC1.rds")




metadata<-dd@meta.data
metadata$Histological<-factor(metadata$Histological,levels=c('metastasis','tumor','normal'),ordered=TRUE)

g<-ggplot(metadata)+
   geom_violin(aes(x=Histological,y=malignant_score,fill=Histological),scale='width')+
   theme_classic()+
   theme(legend.position='none',text=element_text(size=6))+
   scale_fill_manual(values=c('metastasis'='#bb0021','tumor'='#008b45','normal'='#3b4992'))+
   xlab('Histological site')+ylab('malignant score')

ggsave('annotated2/Histological_malignant_score.pdf',plot=g,width=60,height=50,units='mm')



g<-ggplot(metadata)+
   geom_violin(aes(x=Histological,y=non_malignant_score,fill=Histological),scale='width')+
   theme_classic()+
   theme(legend.position='none',text=element_text(size=6))+
   scale_fill_manual(values=c('metastasis'='#bb0021','tumor'='#008b45','normal'='#3b4992'))+
   xlab('Histological site')+ylab('non-malignant score')

ggsave('annotated2/Histological_non_malignant_score.pdf',plot=g,width=60,height=50,units='mm')


metadata<-subset(metadata,subcluster!='Epi_KRT13')
metadata<-subset(metadata,nCount_RNA>3000)

metadata$DE_malignant_log2FC1<-factor(metadata$DE_malignant_log2FC1,levels=c('malignant','non_malignant'),ordered=T)

g<-ggplot(metadata)+
   geom_violin(aes(x=DE_malignant_log2FC1,y=CNV_score,fill=DE_malignant_log2FC1),scale='width')+
   scale_fill_manual(values=c('malignant'='#bb0021','non_malignant'='#3b4992'))+
   theme_classic()+
   theme(legend.position='none',text=element_text(size=6),axis.text=element_text(size=6))+
   xlab('Cell group')+ylab('CNV score')

wilcox.test(metadata$CNV_score[metadata$DE_malignant_log2FC1=='malignant'],metadata$CNV_score[metadata$DE_malignant_log2FC1=='non_malignant'])

ggsave('annotated2/malignant_CNV_score.pdf',plot=g,width=60,height=45,units='mm')

