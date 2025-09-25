
library(ggplot2)
library(Seurat)
library(harmony)

source("Utils.R")


## features
broad_endo<-c("CHGA","CHGB")
gastric_endo<-c("GAST","GHRL","SST")
broad_epi<-c("EPCAM","KRT19","KRT8","KRT7","KRT18","KRT13","KRT6A","KRT5","KRT4","MUC1","CDH1","CLDN4")
gastric_epi<-c("MUC5AC","MUC6","TFF1","TFF2","PGA5","PGA3","PGC","LIPF","ATP4A","ATP4B","CBLIF")
intestinal_epi<-c("REG4","TFF3","SPINK4","MUC2","FABP1","CA1")
stemness<-c("OLFM4","SOX2","LGR5","CCKBR","SOX9","LEFTY1")
fibro<-c("FAP","PDPN","LUM","DCN","COL3A1","COL6A1","COL1A2")
endothelial<-c("PECAM1","VWF","ENG","MCAM","CDH5","PLVAP","ACKR1")
pericyte<-c("RGS5","NOTCH3")
myo<-c("ACTA2","ACTN2","MYL2","MYH2","MYH11","ACTG2")
prol<-c("MKI67","BIRC5","CDK1","HMGB1","STMN1")
hepatocyte<-c("ALB","APOA1","APOA2")
immune<-c("TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2","CD3D","CD3E","CD8A","CD8B",
          "KLRD1","KLRF1",  ## NK?
          "CD79A","CD79B","MS4A1","CD19","IGHD","IGHG1","VPREB3","IGKC","JCHAIN",
		  "XBP1","SDC1","IGHA1","TNFRSF17","MZB1",  ## plasma-related
          "TPSAB1","TPSB2","KIT","CPA3",
		  "CD14","CD163","CD68","FCGR3A","CSF1R","HLA-DQA1","HLA-DPB1","HLA-DRB1","C1QA","C1QB","C1QC","AIF1","TYROBP","FCER1G",
		  "CLEC10A","CD1C","ID2","FCER1A","CLEC4C","LILRA4",  ## DC-related
		  "GZMK","GZMB","GZMM","PRF1","GNLY",
		  "KLRB1","KLRG1","IL7R")
HB<-c("HBA1","HBA2","HBB","HBM","HBD","HBG2","HBG1")

features<-unique(c(broad_endo,gastric_endo,broad_epi,gastric_epi,intestinal_epi,stemness,fibro,endothelial,pericyte,myo,prol,hepatocyte,immune,HB))
## END

d<-readRDS("sep_normalized.rds")

meta<-d@meta.data
meta$Histological[which(meta$Histological=="ptumor")]<-"normal"
meta$Histological[which(meta$Patient %in% c("CR_P1_","CR_P2_","CR_P9_"))]<-"normal"
meta$Histological.site[which(meta$Patient %in% c("CR_P1_","CR_P2_","CR_P9_"))]<-"normal"

d[["Histological"]]<-meta$Histological
d[["Histological.site"]]<-meta$Histological.site

d<-NormalizeData(d)
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


pca.dims<-30

d <- RunTSNE(d, reduction = "harmony", dims = 1:pca.dims)

saveRDS(d,"integrated_harmony.rds")


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


png("./tSNE_results/recluster.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "recluster"))
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
tiff("./tSNE_results/CTmarkerExp_all_logNorm.tif",width = 3800,height = 5000)
print(g)
dev.off()



d <- FindNeighbors(d, reduction = "harmony", dims = 1:pca.dims)
ss<-seq(from=0.1,to=0.3,by=0.1)
for (i in ss){
  print(i)
  d <- FindClusters(d, resolution = i)
  png(paste0("./tSNE_results/tSNE_res_",i,".png"),width = 800,height = 600)
  print(DimPlot(d, reduction = "tsne",label=TRUE))
  dev.off()
}

saveRDS(d,"integrated_harmony.rds")


res<-c(0.1)

for (i in res){
  dir.create(paste0("res_",i))
  print(i)
  d[["seurat_clusters"]]<-d[[paste0("RNA_snn_res.",i)]]
  Idents(d)<-d[["seurat_clusters"]]
  print(table(d[["seurat_clusters"]]))
  print(round(table(d[["seurat_clusters"]])/dim(d)[2]*100,digits=2))
  
## using logCPM normalization data; wilcox method
  celltype.markers<-FindAllMarkers(d,assay="RNA",slot="data",max.cells.per.ident=Inf,only.pos=FALSE)

  celltype.markers<-celltype.markers[,c("cluster","avg_log2FC","pct.1","pct.2","p_val","p_val_adj","gene")]
  celltype.markers<-celltype.markers[which(celltype.markers$p_val_adj<0.05),]

  celltype.markers.neg<-subset(celltype.markers,avg_log2FC<0)
  celltype.markers<-subset(celltype.markers,avg_log2FC>0)

  celltype.markers.neg<-celltype.markers.neg[order(celltype.markers.neg$cluster,celltype.markers.neg$avg_log2FC),]
  write.table(celltype.markers.neg,paste0("res_",i,"/celltype_NEG_markers_Wilcox_all_logNorm.txt"),sep="\t",row.names=FALSE,quote=FALSE)

  celltype.markers<-celltype.markers[order(celltype.markers$cluster,-celltype.markers$avg_log2FC),]
  
  write.table(celltype.markers,paste0("res_",i,"/celltype_markers_Wilcox_all_logNorm.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  celltype.markers.top20<-NULL
  for (j in unique(celltype.markers$cluster)){
    k<-which(celltype.markers$cluster==j)
    if (length(k)>20){
      k<-k[c(1:20)]
    }
    celltype.markers.top20<-rbind(celltype.markers.top20,celltype.markers[k,])
  }
  write.table(celltype.markers.top20,paste0("res_",i,"/celltype_markers_Wilcox_top20_logNorm.txt"),sep="\t",row.names=FALSE,quote=FALSE)

}

d[["main_cluster_res.0.2"]]<-d[["RNA_snn_res.0.2"]]
Idents(d)<-d[["main_cluster_res.0.2"]]

ct<-c("T_NK","Epithelial","Mast","Epithelial","Epithelial",
      "Endocrine","Fibroblast","Epithelial","Proliferative","Epithelial",
	  "Epithelial","Epithelial","Plasma","Epithelial","Myeloid",
	  "B","Fibroblast","Endothelial","Epithelial","T_NK")
names(ct)<-levels(Idents(d))

z<-ct[d[["main_cluster_res.0.2"]]$main_cluster_res.0.2]
d[["main_cluster"]]<-z



saveRDS(d,"pooled_pca.rds")




