source("./Utils.R")

library(ggplot2)


exps<-readRDS('./exp_matrix_gene.rds')
meta<-readRDS('./sample_meta.rds')

exps<-exps[,colnames(exps) %in% meta$GEO_ID]
z<-factor(colnames(exps),levels=meta$GEO_ID,ordered = T)
exps<-exps[,order(z)]


meta$group<-'Low'
meta$group[meta$gene_signature_score>median(meta$gene_signature_score)]<-'High'

meta$CEBPB_exps<-exps['CEBPB',]
saveRDS(meta,'meta_v3.rds')


meta$non_metastasis<-1
meta$non_metastasis[meta$metastasis_liver==1]<-0
meta$non_metastasis[meta$metastasis_other_only==1]<-0

meta$metastasis_group<-"Pri. only"
meta$metastasis_group[meta$metastasis_liver==1]<-"Met. liver"
meta$metastasis_group[meta$metastasis_other_only==1]<-"Met. others"

meta$metastasis_group<-factor(meta$metastasis_group)
g<-ggplot(meta)+
   geom_boxplot(aes(x=metastasis_group,y=CEBPB_exps,fill=metastasis_group),outlier.size = 0.5)+
   theme_classic()+theme(legend.position = 'none',
                         axis.text.y = element_text(size=6),
                         axis.text.x = element_text(size=6,angle=30,vjust=1,hjust=1),
                         axis.title = element_text(size=6))

ggsave(plot=g,"CEBPB_exps_met_groups.pdf",width = 35,height=45,units = 'mm')
