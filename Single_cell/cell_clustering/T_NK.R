
ct<-"T_NK"
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


## marker gene ref:
## 10.1038/s41586-018-0694-x
## 10.1016/j.cell.2017.05.035
## 10.1038/s41467-021-22164-6
## 10.1016/j.jhep.2020.05.039
## 10.1016/j.cell.2020.03.048
## immune-cell-guide
## END

# Hepatocyte<-c("ALB","APOA1","APOC3","ANXA4","DEFB1","AMBP")
Epis<-c("KRT8","KRT18","KRT19","EPCAM")
Fibro<-c("LUM","DCN","COL1A1","COL3A1")
Bc<-c("CD79A","CD79B")
Macro<-c("CD68","CD163","C1QC")
Mono<-c("FCGR3A","FCGR3B","FCN1")

Immune<-c("PTPRC","SRGN")

GdT<-c("TRDC","TRGC1","TRGC2")
Tc<-c("TRAC","TRBC1","TRBC2","CD3D","CD3E")
CD8T<-c("CD8A","CD8B")
CD4T<-c("CD4","IL7R")  ## certain CD8T subclusters, like CD8T-CM, may also express IL7R

Naive<-c("LEF1","SELL","CCR7","TCF7","CD27","CD28","S1PR1")  
EFF<-c("KLRG1","CX3CR1","FCGR3A","PRF1","GZMH","GZMA","GZMB","TBX21","EOMES","GNLY","NKG7","IFNG")
exhausted<-c("LAYN","CXCL13","PDCD1","CTLA4","HAVCR2")
RM<-c("CD69","ITGAE","ZNF683")  ## tissue-resident markers
Circulation<-c("KLF2") ## marker indicating circulation

CD8_RM<-c("CD6","XCL1","XCL2")
CD8_CM<-c("GPR183","CCL5")
CD8_EM<-c("GZMK","CD44","CXCR3","CXCR4")
MAIT<-c("RORA","RORC","ZBTB16","SLC4A10")
IEL<-c("CD160","KIR2DL4","TMIGD2","KLRC1","KLRC2","KLRC3","NR4A1","NR4A2","NR4A3")  ## intraepithelial lymphocytes; also NK cell markers


Treg<-c("FOXP3","IL2RA","IKZF2")  ## Treg and TFH may have close development relationship
TFR<-c("CCR4","IL10")  ## follicular regulatory
Th1<-c("CXCR3","TBX21","IFNG","IL2","LTA")
Th2<-c("GATA3"," PTGDR2","AREG","IL4","IL5","IL10","IL13")
Th17<-c("IL23R","IL17A","CCR6","CCL20","LTB","RORA","KLRB1","FURIN","CTSH","CSF2","MAF","STAT3","IL22")

CD4_CM<-c("ANXA1","ANXA2","PTGER2","ICAM2")
TFH<-c("BCL6","TOX","TOX2","CXCR5","BTLA","ICOS","CD200","PDCD1","ICA1")  ## follicular helper
CD4T_RM<-c("PTGER4","MYADM")

features<-unique(c(Immune,GdT,Tc,CD8T,CD4T,Naive,EFF,exhausted,RM,Circulation,CD8_CM,CD8_EM,MAIT,IEL,TFH,TFR,Th1,Th2,Th17,Treg,CD4_CM,CD4T_RM,CD8_RM,Epis,Fibro,Bc,Macro,Mono))



pca.dims<-30

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


res<-0.9
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
             'CD8T_XCL1','CD4T_CCR7','CD8T_GZMK','CD4T_FOXP3','NK_CD160','CD4T_TNF','CD8T_PLCG2','NK_FCGR3A',
			     'T_NK_MT_ATP6','CD8T_CXCL13','CD8T_KLRB1','CD4T_IL17A','CD4T_BTLA','T_NK_STMN1','CD8T_CCR7','CD4T_MT1X','CD4T_GZMK')
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
