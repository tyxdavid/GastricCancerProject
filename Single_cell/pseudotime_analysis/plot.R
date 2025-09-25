setwd('./main_cluster/subcluster/Epithelial')

library(ggplot2)
library(monocle)
library(Seurat)
library(ggpubr)
library(data.table)
library(ComplexHeatmap)

source('./Utils.R')


## calculate stemness-related score
stemness_markers<-c('EPCAM','ICAM1','THY1','TFRC','CD24','POU5F1','SOX2','PROM1','CD44',
                    'ALDH1A1','BHLHA15','LGR5','AQP5','NANOG','ABCB1','ABCG2','CXCR4','ALCAM','DCLK1','ITGA6')

ST<-readRDS('./main_cluster/subcluster/Epithelial/seurat_for_monocle.rds')

ST_exps<-GetAssayData(ST,assay='RNA')
ST_regulons<-GetAssayData(ST,assay='regulon')

scored_d<-readRDS('./main_cluster/subcluster/Epithelial/Organ_spec_scores_3M_0823/scored_d.rds')

activated_TF<-read.table('./main_cluster/subcluster/Epithelial/scenic/celltype_markers_Wilcox_all.txt',
                         sep="\t",header=TRUE)

TFs<-unlist(read.table('/home/renlab_manager/cisTarget_databases/hg38/allTFs_hg38.txt',sep = "\t",header = F))
organ_genes<-read.table('./GTEx/DEG_GTEx_organs.txt',sep = "\t",header = T)
organ_genes<-organ_genes[organ_genes$avg_log2FC>1 & organ_genes$p_val_adj<0.01 & organ_genes$cluster=='Liver',]
liver_TFs<-intersect(TFs,organ_genes$gene)

liver_TFs<-intersect(liver_TFs,activated_TF)


cytotrace<-readRDS('./main_cluster/subcluster/Epithelial/CytoTRACE_results.rds')
cytotrace<-cytotrace$CytoTRACE

HSMM<-readRDS(paste0('ordered_HSMM.rds'))
genes<-fData(HSMM)$gene_short_name
liver_TFs<-intersect(liver_TFs,genes)
TFs<-intersect(TFs,genes)


g<-plot_cell_trajectory(HSMM, color_by = "subcluster")
ggsave("trajectory.png",plot=g)


g<-plot_cell_trajectory(HSMM, color_by = 'as.factor(State)')
ggsave("cell_state_trajectory.png",plot=g)


g<-plot_cell_trajectory(HSMM, color_by = "Pseudotime")
ggsave("pseudotime.png",plot=g)

g<-plot_cell_trajectory(HSMM, color_by = "Dataset")
ggsave("Dataset.png",plot=g)

g<-plot_cell_trajectory(HSMM, color_by = "DE_malignant_log2FC1")
ggsave("DE_malignant_log2FC1.png",plot=g)

g<-plot_cell_trajectory(HSMM, color_by = "Histological")
ggsave("Histological.png",plot=g)


g<-plot_cell_trajectory(HSMM, color_by = "subcluster") +
    facet_wrap(~subcluster, nrow = 2)
ggsave("trajectory_sep_subclusters.png",plot=g,height=10,width=14)

heatmap_legend_param<-list(title='Expression',grid_height=unit(8,'mm'),grid_width=unit(3,'mm'),title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6))

meta<-pData(HSMM)

meta$differentiation<-cytotrace
meta$stemness<-ST@meta.data$stemness_score1

exps<-as.matrix(HSMM[liver_TFs,])
TF_res<-plot_pseudotime_heatmap_CH(exps=exps,meta=meta,bins=50,cluster_rows=TRUE,show_row_names=TRUE,row_dend_width=unit(3,'mm'),
                                scale_method='none',centering='row',dist_method='euclidean',hclust_method='ward.D2',cluster_num=1,
								heatmap_legend_param=heatmap_legend_param,data_min=-3,data_max=3,smoothed=TRUE,color=rev(rainbow(10)[1:8]))

pdf('liver_TFs.pdf',width=3.937,height=3.15)
print(TF_res$g)
dev.off()


score_mat<-matrix(nrow=4,ncol=nrow(meta))
colnames(score_mat)<-rownames(meta)
rownames(score_mat)<-c('Liver','Stomach','Dedifferentiation','Stemness')
score_mat['Liver',]<-meta$Liver_score1
score_mat['Stomach',]<-meta$Stomach_score1
score_mat['Dedifferentiation',]<-meta$differentiation
score_mat['Stemness',]<-meta$stemness


score_res<-plot_pseudotime_heatmap_CH(exps=score_mat,meta=meta,bins=50,cluster_rows=FALSE,show_row_names=TRUE,
                                scale_method='row',centering='row',
								heatmap_legend_param=heatmap_legend_param,color=rev(rainbow(10)[1:8]),
								smoothed=TRUE)

pdf('Scores.pdf',width=3.937,height=1)
print(score_res$g)
dev.off()




gMet<-ggplot(meta)+
      geom_histogram(aes(x=Pseudotime,fill=Histological),position='fill',bins=50)
gMet<-ggplot_build(gMet)$data[[1]]
gMet<-gMet[,c('fill','count')]
gMet$type<-'tumor'
gMet$type[gMet$fill=='#F8766D']<-'metastasis'

met_mat<-matrix(c(gMet$count[gMet$type=='metastasis'],gMet$count[gMet$type=='tumor']),ncol=2)
rs<-rowSums(met_mat)
for (j in 1:2){
  met_mat[,j]<-met_mat[,j]/rs
}

ha <- HeatmapAnnotation(T_M = anno_barplot(met_mat, gp = gpar(fill = c('#F8766D','#00BFC4'), col = c('#F8766D','#00BFC4'))),height=unit(5,'mm'))

gMag<-ggplot(meta)+
      geom_histogram(aes(x=Pseudotime,fill=DE_malignant_log2FC1),position='fill',bins=50)
gMag<-ggplot_build(gMag)$data[[1]]
gMag<-gMag[,c('fill','count')]
gMag$type<-'non_malignant'
gMag$type[gMag$fill=='#F8766D']<-'malignant'

met_mat<-matrix(c(gMag$count[gMag$type=='malignant'],gMag$count[gMag$type=='non_malignant']),ncol=2)
rs<-rowSums(met_mat)
for (j in 1:2){
  met_mat[,j]<-met_mat[,j]/rs
}

hb <- HeatmapAnnotation(NM = anno_barplot(met_mat, gp = gpar(fill = c('#F8766D','#00BFC4'), col = c('#F8766D','#00BFC4'))),height=unit(5,'mm'))

meta$subcluster_r2<-meta$subcluster
meta$subcluster_r2[! meta$subcluster %in% c('Epi_FABP1','Epi_MDK')]<-'Epi_others'

gCC<-ggplot(meta)+
      geom_histogram(aes(x=Pseudotime,fill=subcluster_r2),position='fill',bins=50)

# ggsave(plot=gCC,'gCC_preview.png')

gCC<-ggplot_build(gCC)$data[[1]]
gCC<-gCC[,c('fill','count')]
gCC$type<-'Epi_others'
gCC$type[gCC$fill=='#F8766D']<-'Epi_FABP1'
gCC$type[gCC$fill=='#00BFC4']<-'Epi_MDK'

met_mat<-matrix(c(gCC$count[gCC$type=='Epi_FABP1'],gCC$count[gCC$type=='Epi_MDK'],gCC$count[gCC$type=='Epi_others']),ncol=3)
rs<-rowSums(met_mat)
for (j in 1:3){
  met_mat[,j]<-met_mat[,j]/rs
}

hc <- HeatmapAnnotation(CC = anno_barplot(met_mat, gp = gpar(fill = c('#F8766D','#00BFC4','grey50'), col = c('#F8766D','#00BFC4','grey50'))),height=unit(5,'mm'))


ht_list <- TF_res$g %v% score_res$g %v% ha %v% hb %v% hc

pdf('liver_TFs_score_merged_1025.pdf',width=3.5433,height=3.1496)
draw(ht_list)
dev.off()


ggtheme<-theme(axis.text = element_text(size=6),
               text=element_text(size=6),axis.title=element_text(size=6),
               legend.text=element_text(size=6),legend.title=element_text(size=6),
               legend.margin=margin(0,0,0,0))


g1<-plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size=0.1,show_branch_points=FALSE,cell_link_size=0.5)+ggtheme+theme(legend.position='top',legend.key.width=unit(4,'mm'),,legend.key.height=unit(2,'mm'))
ggsave("trajectory_pseudotime.pdf",plot=g1,width=45,height=50,units='mm')


gts<-plot_cell_trajectory(HSMM, color_by = "subcluster",cell_size=0.01,cell_link_size=0.3,show_branch_points=FALSE) +
    facet_wrap(~subcluster, nrow = 2)+ggtheme+theme(legend.position='none',axis.ticks=element_blank(),axis.text=element_blank())
ggsave("trajectory_sep_subclusters.pdf",plot=gts,height=50,width=70,units='mm')



# ## distance along the trajectory

dist_pesudotime<-dist(as.matrix(meta$Pseudotime,ncol=1),diag=TRUE)
dist_pesudotime<-as.matrix(dist_pesudotime)

clusters<-unique(meta$subcluster)

res<-NULL
x1<-which((meta$subcluster == 'Epi_FABP1'))

for (i in clusters){
  x2<-which((meta$subcluster == i))
  res_sub<-data.frame(C2=rep(i,length(x2)),dis=as.numeric(dist_pesudotime[x2,x1]))
  res<-rbind(res,res_sub)
}

res<-subset(res,C2!='Epi_FABP1')
res$C2<-factor(res$C2,levels=c('Epi_MDK','Epi_CD74','Epi_HMGB2','Epi_SOX4','Epi_MUC2'),ordered=TRUE)

ggtheme<-theme(text  = element_text(size=6),
               title  = element_text(size=6),
               axis.text.x=element_text(size=6,angle=45,vjust=1,hjust=1))

g<-ggplot(res,aes(x=C2,y=dis,fill=C2))+
  geom_boxplot(outlier.size = 0.1,width=0.35)+
  theme_classic()+ggtheme+theme(legend.position = 'none')+
  ylab('Distance to cells of Epi_FABP1 cluster')+xlab('Cell cluster')+
  scale_fill_hue()

ggsave('Distance_to_FABP1.pdf',width=47,height=51,units='mm')

aov_res<-aov(dis ~ C2, data = res)
print(summary(aov_res))
print(TukeyHSD(aov_res))


dist_matrix<-matrix(ncol=length(clusters),nrow=length(clusters))
dimnames(dist_matrix)<-list(clusters,clusters)

for (i in 1:length(clusters)){
  for (j in i:length(clusters)){
    x1<-which((meta$subcluster %in% clusters[j]))
	x2<-which((meta$subcluster %in% clusters[i]))
	dis<-sum(dist_pesudotime[x2,x1])/(length(x1)*length(x2))
    dist_matrix[i,j]<-dis
  }
}

saveRDS(dist_matrix,'dist_matrix.rds')



dist_matrix_melt<-reshape2::melt(dist_matrix)

dist_matrix_melt$value<-round(dist_matrix_melt$value,1)

g<-ggplot(dist_matrix_melt)+
   geom_tile(aes(x=Var1,y=Var2,fill=value))+
   scale_fill_gradient2(low='red',high='blue',mid='#EEEEEE',midpoint=mean(dist_matrix,na.rm=TRUE),name='Distance',na.value='grey80')+
   geom_text(aes(x=Var1,y=Var2,label=value),color='black',size=2)+
   xlab('Cluster 1')+ylab('Cluster 2')+
   theme_classic()+
   theme(text=element_text(size=6),axis.text.y=element_text(size=6),axis.text.x=element_text(size=6,angle=45,vjust=1,hjust=1),axis.title=element_text(size=6),
           legend.text=element_text(size=6),legend.title=element_text(size=6),
			   legend.margin=margin(0,0,0,0),legend.key.width=unit(3,'mm'),legend.key.height=unit(3,'mm'))

ggsave("pseudotime_dist_subclusters.pdf",plot=g,height=50,width=60,units='mm')


# ## liver score per cluster

meta$subcluster<-factor(meta$subcluster,levels=c('Epi_MDK','Epi_FABP1','Epi_SOX4','Epi_MUC2','Epi_CD74','Epi_HMGB2'),ordered=TRUE)

print(wilcox.test(meta$Liver_score1[meta$subcluster=='Epi_MDK'],meta$Liver_score1[meta$subcluster=='Epi_FABP1']))
print(wilcox.test(meta$Liver_score1[meta$subcluster=='Epi_MDK'],meta$Liver_score1[meta$subcluster=='Epi_SOX4']))
print(wilcox.test(meta$Liver_score1[meta$subcluster=='Epi_FABP1'],meta$Liver_score1[meta$subcluster=='Epi_SOX4']))


g<-ggplot(meta)+
   geom_boxplot(aes(x=subcluster,y=Liver1,fill=subcluster),outlier.size = 0.5,outlier.stroke = 0)+
   scale_fill_hue()+
   theme_classic()+ggtheme+theme(legend.position='none',axis.text.x=element_text(size=6,angle=45,vjust=1,hjust=1))+
   ylab('Liver module expression')+xlab('Cell cluster')

ggsave(plot=g,'liver_per_cluster.pdf',height=50,width=50,units='mm')



## CEBPB exps per cluster
ST<-readRDS('./main_cluster/subcluster/Epithelial/seurat_for_monocle.rds')
exps<-GetAssayData(ST)['CEBPB',]
dc<-data.frame(cluster=meta$subcluster,CEBPB=exps)


print(wilcox.test(dc$CEBPB[dc$cluster=='Epi_MDK'],dc$CEBPB[dc$cluster=='Epi_FABP1']))
print(wilcox.test(dc$CEBPB[dc$cluster=='Epi_MDK'],dc$CEBPB[dc$cluster=='Epi_SOX4']))
print(wilcox.test(dc$CEBPB[dc$cluster=='Epi_FABP1'],dc$CEBPB[dc$cluster=='Epi_SOX4']))



## plot SCENIC activate TFs heatmap

ST_regulons<-GetAssayData(ST,assay='regulon')

meta<-ST@meta.data

activated_TF<-read.table('./main_cluster/subcluster/Epithelial/scenic/celltype_markers_Wilcox_all.txt',
                         sep="\t",header=TRUE)
activated_TF<-subset(activated_TF,avg_log2FC>0.01 & p_val_adj<0.01)
cluster_ord<-unique(activated_TF$cluster)
activated_TF<-unique(activated_TF$gene)

organ_genes<-read.table('./GTEx/DEG_GTEx_organs.txt',sep = "\t",header = T)
organ_genes<-organ_genes[organ_genes$avg_log2FC>1 & organ_genes$p_val_adj<0.01 & organ_genes$cluster=='Liver',]
liver_TFs<-intersect(activated_TF,organ_genes$gene)


ST_regulons<-ST_regulons[activated_TF,]

res<-plot_heatmap_avg_exps(as.matrix(ST_regulons),class_label = as.character(meta$subcluster),
                           scaling = T,centering = T,
                           high_col='#ff0e04',low_col='#11098c',mid_col='#ffffff',
                           cluster_order = cluster_ord)
g<-res$g
g<-g+theme(axis.ticks=element_line(linewidth=0.25),
        axis.text.y=element_text(size=6),
        axis.text.x.top=element_text(size=6,angle=45,hjust=0,vjust=0),
		   axis.title=element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6),
		   legend.key.height=unit(5,'mm'),legend.key.width=unit(2,'mm'))+
		   xlab('Cell cluster')+ylab('TF regulons')

show_label_func<-function(x,included_labels=NULL){
  output<-x
  z<-which(! output %in% included_labels)
  output[z]<-rep('',length(z))
  names(output)<-x
  return(output)
}

z<-show_label_func(rownames(res$avg_exps),included_labels=liver_TFs)

g<-g+scale_y_discrete(labels=z)

ggsave(plot=g,'regulons_HM_v2.pdf',units='mm',width=60,height=75)
