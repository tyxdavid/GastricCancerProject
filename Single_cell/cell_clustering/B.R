
ct<-"B"

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
## 10.1016/j.it.2022.01.003
## 10.1126/sciimmunol.abe6291
## 10.1016/j.celrep.2022.110713
## 10.1016/j.cell.2020.03.048
## immune-cell-guide
## END

Epis<-c("KRT8","KRT18","KRT19","EPCAM")
Bc<-c("CD79A","CD79B")
Tc<-c("CD3D","CD3E","TRAC")

Immune<-c("PTPRC","SRGN")
Antigen_process<-c("HLA-DQA1","CD74")

IGs<-c("IGHM","IGHD","IGHA1","IGHG1","IGKC","JCHAIN")
preB<-c("IL7R","LEF1","SOX4","PAX5","RAG1","RAG2","IGLL1")
immatureB<-c("CD19","CD24","EBF1","EIF4EBP1","VPREB3")
memoryB<-c("FCRL4","CCR1","GPR183","FCER2","CD44","KLF2")
GCB<-c("AICDA","GCSAM","RGS13","BCL6")
plasma<-c("SDC1","MZB1","ATF4")

features<-unique(c(Immune,Antigen_process,Bc,IGs,preB,immatureB,memoryB,GCB,plasma,Epis,Tc))



pca.dims<-15

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

g<-FeaturePlot(d, features = features,reduction = "tsne",label = FALSE,ncol=7,slot="data")
png("./tSNE_results/CTmarkerExp_all_logNorm.png",width = 3000,height = 3000)
print(g)
dev.off()



d <- FindNeighbors(d, reduction = "harmony", dims = 1:pca.dims)
ss<-seq(from=0.1,to=1,by=0.1)
for (i in ss){
  print(i)
  d <- FindClusters(d, resolution = i)

  png(paste0("./tSNE_results/tSNE_res_",i,".png"),width = 800,height = 600)
  print(DimPlot(d, reduction = "tsne",label=TRUE))
  dev.off()
}


saveRDS(d,"integrated_harmony.rds")


res<-0.4
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
ggsave("DotPlot_markers.png",width=10,height=6,plot=g)



d <- RenameIdents(d, 'B_GPR183','B_IGHD','B_EGR1','B_MKI67','B_MEF2B',
                  'B_NME1','B_MT_CO2','B_FCRL4')
## END



d[["subcluster"]]<-Idents(d)



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

metadata<-d@meta.data
enriched_idents(metadata,idents_tested="Histological",idents_based="subcluster",tmpdir="enriched_idents_Histological_subcluster",angle=45)
enriched_idents(metadata,idents_tested="subcluster",idents_based="Histological",tmpdir="enriched_idents_subcluster_Histological",angle=45)


saveRDS(d,"integrated_harmony.rds")



features<-c('PTPRC','SRGN',
            'CD79A','CD79B',
            "CCR1","GPR183","FCER2","CD44","KLF2",
			'CD70',
			'NME1','NME2','ENO1','FABP5','MIR155HG',
			'FCRL4','DUSP4','FCRL5','ITGAX',
			'IGHD','IGHM','YBX3','KLF3',
			"AICDA","GCSAM","RGS13","BCL6",'MEF2B','LRMP','HOPX',
			'EGR1','CD69','DUSP2','JUN')

proliferative<-c("TYMS","MKI67","STMN1")


features<-unique(c(features,proliferative))

g<-DotPlot(d, 
      features = features, 
      cols = c('lightgrey', 'red')) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
ggsave("DotPlot_markers_annotated.png",width=18,height=6,plot=g)

