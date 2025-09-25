setwd('./GTEx')
library(parallel)
library(GSVA)
library(ggplot2)
source('./Utils.R')

tpm<-readRDS('./GTEx/exps_log2TPM_all_samples.rds')
meta<-readRDS('./GTEx/meta_all_samples.rds')

DEGs<-readRDS("./GTEx/wilcox_all_organs/DEG_GTEx_organs.rds")
DEGs<-subset(DEGs, avg_log2FC>1 & p_val_adj<0.01)

print(table(DEGs$cluster))

organ_per_genes<-table(DEGs$gene)
retain_genes<-as.integer(organ_per_genes)
names(retain_genes)<-names(organ_per_genes)
retain_genes<-retain_genes[retain_genes<=3]

DEGs<-DEGs[DEGs$gene %in% names(retain_genes),]


write.table(DEGs[order(DEGs$cluster,-DEGs$avg_log2FC),],'selected_DEGs.txt',row.names=F,quote=F,sep="\t")

print(table(DEGs$cluster))

geneSets<-list()

for (i in unique(DEGs$cluster)){
  geneSets[[gsub(" ","_",i,fixed=TRUE)]]<-DEGs$gene[which(DEGs$cluster==i)]
}


set.seed(1775)

ssgsea<-gsva(tpm,geneSets,method="ssgsea",parallel.sz=64)

ssgsea<-t(ssgsea)
meta<-cbind(meta,ssgsea)

saveRDS(ssgsea,"GTEx_organ_scores.rds")

colnames(meta)[(ncol(meta)-ncol(ssgsea)+1):ncol(meta)]<-paste0("Score_",colnames(meta)[(ncol(meta)-ncol(ssgsea)+1):ncol(meta)])

saveRDS(meta,'meta_all_samples.rds')

ssgsea<-t(ssgsea)

res<-plot_heatmap_avg_exps(exps=ssgsea,class_label=meta$SMTS,scaling=TRUE,ignore_class_sample_size=FALSE)

res$g<-res$g+theme(aspect.ratio=1,axis.text.x.top=element_text(angle=90,vjust=0.5,size=6,hjust=0),axis.text.y=element_text(size=6),axis.title=element_text(size=6),
                   legend.key.width=unit(3,'mm'),legend.key.height=unit(6,'mm'),legend.text=element_text(size=6),legend.title=element_text(size=6),legend.position='right',
				   axis.ticks=element_line(linewidth=0.2))+
	   scale_fill_gradient(low='white',high='red')+
	   xlab('sampled organs')+ylab('expression profiles (defined by organ-specific DEGs)')
ggsave(plot=res$g,'GTEx_organ_expression_score.pdf',height=90,width=100,units='mm')


