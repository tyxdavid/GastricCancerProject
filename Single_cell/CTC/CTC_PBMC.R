
library(Seurat)
library(ggplot2)
library(DESeq2)

source('Utils.R')

##################
##### Fig 1e #####
##################


exps<-readRDS('aggr_CTC_PBMC_count.rds')

gL<-read.table('geneLength.txt',header=TRUE,sep="\t")
tpm<-log1p(TPM(exps,gL))

meta<-data.frame(Histological=c(rep('CTC',27),rep('PBMC',13)),
                 Histological.site=c(rep('CTC',27),rep('PBMC',13)),
				 Dataset=c(rep('US',27),rep('GSE107011',13)),
				 Patient=substr(colnames(exps),1,6))

rownames(meta)<-colnames(exps)

d<-CreateSeuratObject(counts=exps,meta.data=meta)

d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")

d<-SetAssayData(d, slot = "data", new.data=tpm)

d <- FindVariableFeatures(d, selection.method = "vst", nfeatures = 2000)
d <- ScaleData(d)
d <- RunPCA(d, features = VariableFeatures(object = d),npcs=30)

png("elbowplot_pca.png",width=800,height=600)
print(ElbowPlot(d,ndims=30,reduction="pca"))
dev.off()

pca.dims<-5
d <- RunTSNE(d, reduction = "pca", dims = 1:pca.dims,reduction.name="tsne",perplexity = 5)



dir.create("tSNE_PCA_results")

png("./tSNE_PCA_results/Patient.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "Patient"))
dev.off()

png("./tSNE_PCA_results/Histological.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "Histological"))
dev.off()


png("./tSNE_PCA_results/Dataset.png",width = 800,height = 600)
print(DimPlot(d, reduction = "tsne",group.by = "Dataset"))
dev.off()


g<-FeaturePlot(d, features = "nCount_RNA",reduction = "tsne",label = FALSE)
png("./tSNE_PCA_results/RNA_counts.png",width = 800,height = 600)
print(g)
dev.off()

g<-FeaturePlot(d, features = "nFeature_RNA",reduction = "tsne",label = FALSE)
png("./tSNE_PCA_results/Features_counts.png",width = 800,height = 600)
print(g)
dev.off()

g<-FeaturePlot(d, features = "percent.mt",reduction = "tsne",label = FALSE)
png("./tSNE_PCA_results/percent_mito.png",width = 800,height = 600)
print(g)
dev.off()



d <- FindNeighbors(d, reduction = "pca", dims = 1:pca.dims)
ss<-seq(from=0.1,to=1.5,by=0.1)
for (i in ss){
  print(i)
  d <- FindClusters(d, resolution = i)
 
  png(paste0("./tSNE_PCA_results/tSNE_res_",i,".png"),width = 800,height = 600)
  print(DimPlot(d, reduction = "tsne",label=TRUE))
  dev.off()
}



res<-1.2
i<-res
dir.create(paste0("res_",i))
  print(i)
  d[["seurat_clusters"]]<-d[[paste0("RNA_snn_res.",i)]]
  Idents(d)<-d[["seurat_clusters"]]
  print(table(d[["seurat_clusters"]]))
  print(round(table(d[["seurat_clusters"]])/dim(d)[2]*100,digits=2))
  
  saveRDS(d,"d.rds")

  
  celltype.markers<-FindAllMarkers(d,assay="RNA",max.cells.per.ident=Inf,only.pos=FALSE,test.use='DESeq2')

  celltype.markers<-celltype.markers[,c("cluster","avg_log2FC","pct.1","pct.2","p_val","p_val_adj","gene")]
  celltype.markers<-celltype.markers[which(celltype.markers$p_val_adj<0.05),]

  celltype.markers.neg<-subset(celltype.markers,avg_log2FC<0)
  celltype.markers<-subset(celltype.markers,avg_log2FC>0)

  celltype.markers.neg<-celltype.markers.neg[order(celltype.markers.neg$cluster,celltype.markers.neg$avg_log2FC),]
  write.table(celltype.markers.neg,paste0("res_",i,"/celltype_NEG_markers_DESeq2_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)

  celltype.markers<-celltype.markers[order(celltype.markers$cluster,-celltype.markers$avg_log2FC),]
  
  write.table(celltype.markers,paste0("res_",i,"/celltype_markers_DESeq2_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  celltype.markers.top20<-NULL
  for (j in unique(celltype.markers$cluster)){
    k<-which(celltype.markers$cluster==j)
    if (length(k)>20){
      k<-k[c(1:20)]
    }
    celltype.markers.top20<-rbind(celltype.markers.top20,celltype.markers[k,])
  }
  write.table(celltype.markers.top20,paste0("res_",i,"/celltype_markers_DESeq2_top20.txt"),sep="\t",row.names=FALSE,quote=FALSE)


g<-DimPlot(d, reduction = "tsne",label=FALSE,pt.size=1.5,label.size=3)
ggsave(plot=g,'annotated/tSNE_annotated.pdf',height=7,width=8)
