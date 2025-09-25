
ct<-"Epithelial"
dir.create(ct)
setwd(ct)

library(ggplot2)
library(Seurat)
source("Utils.R")

library(harmony)
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


features = unique(c("EPCAM","KRT8","KRT18","KRT19","MUC1",  ## Pan-epithelial
	               "TFF1","TFF2","GKN1","GKN2","MUC5AC","MUCL3", 
				   "MUC6","PRR4","LTF","LYZ",  
				   "PGA3","PGA5","PGA4","LIPF","PGC","MT1F","MT1M",  
				   "CKB","ATP4B","ATP4A","CBLIF",  
				   "FABP1","ALDOB","AGPAT2","CDHR5","APOA1","APOC3","PCK1","TM4SF4","MTTP", 
				   "SPINK4","MUC2","TFF3","ZG16","ITLN1","FCGBP", 
				   "H2AFZ","HMGB2","CENPF","MKI67",  
				   "KRT13","SPRR3","KRT5","CRNN","CRCT1","SPINK5","CSTB","S100A9","S100A8" 
				     ))


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

g<-FeaturePlot(d, features = features,reduction = "tsne",label = FALSE,ncol=9,slot="data")
png("./tSNE_results/CTmarkerExp_all_logNorm.png",width = 3800,height = 5000)
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


res<-0.3
i<-res
dir.create(paste0("res_",i))
  print(i)
  d[["seurat_clusters"]]<-d[[paste0("RNA_snn_res.",i)]]
  Idents(d)<-d[["seurat_clusters"]]
  print(table(d[["seurat_clusters"]]))
  print(round(table(d[["seurat_clusters"]])/dim(d)[2]*100,digits=2))
  
# saveRDS(d,"integrated_harmony.rds")

  
## using logCPM normalization data; wilcox method
  celltype.markers<-FindAllMarkers(d,assay="RNA",slot="data",max.cells.per.ident=Inf,only.pos=FALSE)

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
ggsave("DotPlot_markers.png",width=18,height=6,plot=g)



d <- RenameIdents(d,
             'Epi_MUC5AC','Epi_SOX4','Epi_HMGB2','Epi_CD24','Epi_MUC6',
			    'Epi_FABP1','Epi_MT_ND2','Epi_SPARC','Epi_KRT13',
			     'Epi_PGA5','Epi_CD74','Epi_MUC2','Epi_CKB','Epi_MDK')
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
