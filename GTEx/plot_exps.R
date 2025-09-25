setwd('./GTEx')
library(data.table)
library(ggplot2)
source('./Utils.R')

exps_meta<-fread('rna_tissue_gtex.tsv',header = T,sep = "\t")


liver_TFs<-c('CEBPB','BCL3','HNF4A','ONECUT2','ATF5','IRF7','CEBPA','NR1H3','NR1I2','NR1H4','AR')

exps_meta<-exps_meta[`Gene name` %in% liver_TFs,]

exps_meta$Tissue<-NULL
exps_meta$TPM<-NULL
exps_meta$nTPM<-NULL
exps_meta$pTPM<-NULL
exps_meta$Sample_counts<-NULL
exps_meta<-unique(exps_meta)

exps_meta$TPM_new<-log2(exps_meta$TPM_new+1e-2)

exps_meta$Gene<-exps_meta$`Gene name`
exps_meta$`Gene name`<-NULL

exps_meta$exps<-exps_meta$TPM_new
exps_plot<-exps_meta[,c('Gene','Tissue','exps')]

exps_plot<-dcast(exps_plot,Tissue~Gene,value.var = 'exps')
exps_plot<-as.matrix(exps_plot,rownames = 'Tissue')
exps_plot<-scale(exps_plot,scale = F)
exps_plot<-apply(exps_plot,2,norm01)

exps_plot<-melt(exps_plot)

sort_value<-exps_plot[exps_plot$Var2=='CEBPB',]
exps_plot$Var1<-factor(exps_plot$Var1,levels=sort_value$Var1[order(sort_value$value,decreasing = T)],ordered = T)
exps_plot<-exps_plot[order(exps_plot$Var1),]
exps_plot$Var2<-factor(exps_plot$Var2,levels=rev(liver_TFs),ordered = T)

md<-(max(exps_plot$value)-min(exps_plot$value))/2+min(exps_plot$value)
g<-ggplot(exps_plot)+
    geom_tile(aes(x=Var1,y=Var2,fill=value))+
    scale_fill_gradient(high='red',low='white')+
    theme_classic()+
    scale_x_discrete(position = "top")+
    theme(axis.text.x.top=element_text(angle=90,vjust=0.5,hjust=0,size=6),legend.position = 'none',
          axis.text.y=element_text(size=6))+
   xlab('Sampled organs')+ylab('TFs')

ggsave('liver_11TFs_exps.pdf',width=76,height=58,units='mm')
