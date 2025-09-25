
ct<-"Myeloid"

dir.create(ct)
setwd(ct)

library(ggplot2)
library(Seurat)

library(harmony)

source("Utils.R")


d<-readRDS("./main_cluster/pooled_pca.rds")
d<-subset(d,main_cluster==ct)

d[["umap_pca_main_cluster"]]<-d[["umap_pca"]]
d[["tsne_pca_main_cluster"]]<-d[["tsne_pca"]]

d[["umap_pca"]]<-NULL
d[["tsne_pca"]]<-NULL
d[["umap"]]<-NULL
d[["tsne"]]<-NULL

gc()

d<-FindVariableFeatures(d)
d<-ScaleData(d,vars.to.regress=c("percent.mt","nCount_RNA"))

d <- RunPCA(d, features = VariableFeatures(object = d))
d <- RunHarmony(d, c("Dataset","Patient"),max.iter.harmony=50)


png("elbowplot_pca.png",width=800,height=600)
print(ElbowPlot(d,ndims=50,reduction="pca"))
dev.off()

png("elbowplot_harmony.png",width=800,height=600)
print(ElbowPlot(d,ndims=50,reduction="harmony"))
dev.off()

saveRDS(d,"integrated_harmony.rds")




## marker gene ref:
## 10.1016/j.cell.2021.01.010
## 10.1016/j.immuni.2019.03.009
## immune-cell-guide
## END

Epis<-c("KRT8","KRT18","KRT19","EPCAM","MUC5AC")
Bc<-c("CD79A","CD79B","MS4A1","MZB1","IGHA1","IGHM")
Tc<-c("CD3D","CD3E","TRAC")
Fibro<-c("TAGLN","LUM","DCN","COL3A1")

Immune<-c("PTPRC","SRGN")
Antigen_process<-c("HLA-DQA1","CD74")

pDCs<-c("LILRA4","GZMB","IL3RA")
cDCs<-c("CLEC9A","FLT3","IDO1","CD1C","FCER1A","LAMP3","CCR7","FSCN1")
Monos<-c("CD14","FCN1","S100A9","S100A8","FCGR3A","LST1","LILRB2")
Mcs<-c("CD68","CD163","MRC1","C1QC","C1QA","GPNMB","APOE","APOC1","SPP1","LYVE1","SELENOP","INHBA","IL1RN","CCL4","NLRP3","EREG","IL1B","ISG15","FN1")

features<-unique(c(Immune,Antigen_process,pDCs,cDCs,Monos,Mcs,Epis,Fibro,Bc,Tc))




dir.create("tSNE_PCA_results")

png("./tSNE_PCA_results/Patient.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne_pca",group.by = "Patient"))
dev.off()

png("./tSNE_PCA_results/Histological.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne_pca",group.by = "Histological"))
dev.off()

png("./tSNE_PCA_results/Histological_site.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne_pca",group.by = "Histological.site"))
dev.off()


png("./tSNE_PCA_results/Dataset.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne_pca",group.by = "Dataset"))
dev.off()


g<-FeaturePlot(d, features = "nCount_RNA",reduction = "tsne_pca",label = FALSE)
png("./tSNE_PCA_results/UMI_counts.png",width = 800,height = 600)
print(g)
dev.off()

g<-FeaturePlot(d, features = "nFeature_RNA",reduction = "tsne_pca",label = FALSE)
png("./tSNE_PCA_results/Features_counts.png",width = 800,height = 600)
print(g)
dev.off()

g<-FeaturePlot(d, features = "percent.mt",reduction = "tsne_pca",label = FALSE)
png("./tSNE_PCA_results/percent_mito.png",width = 800,height = 600)
print(g)
dev.off()


g<-FeaturePlot(d, features = features,reduction = "tsne_pca",label = FALSE,ncol=8,slot="data")
png("./tSNE_PCA_results/CTmarkerExp_all_logNorm.png",width = 3000,height = 3000)
print(g)
dev.off()


pca.dims<-20

d <- RunTSNE(d, reduction = "harmony", dims = 1:pca.dims)


dir.create("tSNE_results")

png("./tSNE_results/Patient.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "Patient"))
dev.off()

png("./tSNE_results/Histological.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "Histological"))
dev.off()

png("./tSNE_results/Histological_site.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "Histological.site"))
dev.off()


png("./tSNE_results/Dataset.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "Dataset"))
dev.off()

g<-FeaturePlot(d, features = "nCount_RNA",reduction = "tsne",label = FALSE)
png("./tSNE_results/UMI_counts.png",width = 800,height = 600)
print(g)
dev.off()

g<-FeaturePlot(d, features = "nFeature_RNA",reduction = "tsne",label = FALSE)
png("./tSNE_results/Features_counts.png",width = 800,height = 600)
print(g)
dev.off()

g<-FeaturePlot(d, features = "percent.mt",reduction = "tsne",label = FALSE)
png("./tSNE_results/percent_mito.png",width = 800,height = 600)
print(g)
dev.off()

g<-FeaturePlot(d, features = features,reduction = "tsne",label = FALSE,ncol=8,slot="data")
png("./tSNE_results/CTmarkerExp_all_logNorm.png",width = 3000,height = 3000)
print(g)
dev.off()



d <- FindNeighbors(d, reduction = "harmony", dims = 1:pca.dims)
ss<-seq(from=0.1,to=1.5,by=0.1)
for (i in ss){
  print(i)
  d <- FindClusters(d, resolution = i)

  png(paste0("./tSNE_results/tSNE_res_",i,".png"),width = 800,height = 600)
  print(DimPlot(d, reduction = "tsne",label=TRUE))
  dev.off()
}


saveRDS(d,"integrated_harmony.rds")


res<-0.8
i<-res
dir.create(paste0("res_",i))
  print(i)
  d[["seurat_clusters"]]<-d[[paste0("RNA_snn_res.",i)]]
  Idents(d)<-d[["seurat_clusters"]]
  print(table(d[["seurat_clusters"]]))
  print(round(table(d[["seurat_clusters"]])/dim(d)[2]*100,digits=2))
  

## using normalization data
  celltype.markers<-FindAllMarkers(d,assay="RNA",slot="data",max.cells.per.ident=Inf,only.pos=FALSE,test.use='MAST')

  celltype.markers<-celltype.markers[,c("cluster","avg_log2FC","pct.1","pct.2","p_val","p_val_adj","gene")]
  celltype.markers<-celltype.markers[which(celltype.markers$p_val_adj<0.05),]

  celltype.markers.neg<-subset(celltype.markers,avg_log2FC<0)
  celltype.markers<-subset(celltype.markers,avg_log2FC>0)

  celltype.markers.neg<-celltype.markers.neg[order(celltype.markers.neg$cluster,celltype.markers.neg$avg_log2FC),]
  write.table(celltype.markers.neg,paste0("res_",i,"/celltype_NEG_markers_Wilcox_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)

  celltype.markers<-celltype.markers[order(celltype.markers$cluster,-celltype.markers$avg_log2FC),]
  
  write.table(celltype.markers,paste0("res_",i,"/celltype_markers_Wilcox_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  celltype.markers.top20<-NULL
  for (j in unique(celltype.markers$cluster)){
    k<-which(celltype.markers$cluster==j)
    if (length(k)>20){
      k<-k[c(1:20)]
    }
    celltype.markers.top20<-rbind(celltype.markers.top20,celltype.markers[k,])
  }
  write.table(celltype.markers.top20,paste0("res_",i,"/celltype_markers_Wilcox_top20.txt"),sep="\t",row.names=FALSE,quote=FALSE)



g<-DotPlot(d, 
      features = features, 
      cols = c('lightgrey', 'red')) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
ggsave("DotPlot_markers.png",width=10,height=6,plot=g)

metadata<-d@meta.data
enriched_idents(metadata,idents_tested="Histological",idents_based="seurat_clusters",tmpdir="enriched_idents_Histological_cluster",angle=90)
enriched_idents(metadata,idents_tested="seurat_clusters",idents_based="Histological",tmpdir="enriched_idents_cluster_Histological",angle=90)






d <- RenameIdents(d, 'Mc_C1QC','Mono','Mc_APOE','cDC_CD1C','Mc_INHBA','Mc_LYVE1','Mc_FN1','Myeloid_MKI67','Mc_ISG15','cDC_LAMP3',
                  'Neutro','cDC_CLEC9A','pDC_LILRA4')
## END


d[["subcluster"]]<-Idents(d)

saveRDS(d,"integrated_harmony.rds")


dir.create("annotated")
  
## using logCPM normalization data; wilcox method
  celltype.markers<-FindAllMarkers(d,assay="RNA",slot="data",max.cells.per.ident=Inf,only.pos=FALSE)

  celltype.markers<-celltype.markers[,c("cluster","avg_log2FC","pct.1","pct.2","p_val","p_val_adj","gene")]
  celltype.markers<-celltype.markers[which(celltype.markers$p_val_adj<0.05),]

  celltype.markers.neg<-subset(celltype.markers,avg_log2FC<0)
  celltype.markers<-subset(celltype.markers,avg_log2FC>0)

  celltype.markers.neg<-celltype.markers.neg[order(celltype.markers.neg$cluster,celltype.markers.neg$avg_log2FC),]
  write.table(celltype.markers.neg,paste0("annotated","/celltype_NEG_markers_Wilcox_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)

  celltype.markers<-celltype.markers[order(celltype.markers$cluster,-celltype.markers$avg_log2FC),]
  
  write.table(celltype.markers,paste0("annotated","/celltype_markers_Wilcox_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  celltype.markers.top20<-NULL
  for (j in unique(celltype.markers$cluster)){
    k<-which(celltype.markers$cluster==j)
    if (length(k)>20){
      k<-k[c(1:20)]
    }
    celltype.markers.top20<-rbind(celltype.markers.top20,celltype.markers[k,])
  }
  write.table(celltype.markers.top20,paste0("annotated","/celltype_markers_Wilcox_top20.txt"),sep="\t",row.names=FALSE,quote=FALSE)


png(paste0("./annotated/tSNE_annotated.png"),width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",label=TRUE))
dev.off()



features<-c('PTPRC','SRGN',
            'CD68','CD163','MRC1','AIF1',
			'C1QC','C1QA','C1QB','SELENOP','SLC40A1',
			'APOE','APOC1','GPNMB','CCL18','MMP12','MMP9',
			'LYVE1','RNASE1','MARCO','PLTP','KLF2',
			'ISG15','CXCL10','CXCL9','CCL8','IFIT2',
			'INHBA','IL1RN','IL1B','CCL4','CCL20','CXCL5','CXCL3','CXCL1',
			'FN1','FCGR2B','PLTP',
			'CD14','FCN1','S100A8','S100A9',
			'FCGR3A','FCGR3B','LST1','LILRB2',
			'HLA-DQA1','CD74',
			'CD1C','FCER1A','CLEC10A',
			'CLEC9A','FLT3','IDO1',
			'LAMP3','CCR7','FSCN1',
			'LILRA4','GZMB','IL3RA')

proliferative<-c("TYMS","MKI67","STMN1")


features<-unique(c(features,proliferative))

g<-DotPlot(d, 
      features = features, 
      cols = c('lightgrey', 'red')) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
ggsave("DotPlot_markers_annotated.png",width=18,height=6,plot=g)



