setwd("./main_cluster/subcluster")

ct<-"Endocrine"

dir.create(ct)
setwd(ct)

library(ggplot2)
library(Seurat)

source("./Utils.R")

d<-subset(d,seurat_clusters %in% c(0,1,3:6))


d <- RenameIdents(d, '0' = 'G_cell',
                     '1' = '1',
					 '3' = 'D_cell',
					 '4' = 'EC_cell_like',
					 '5' = 'EC_cell',
					 '6' = 'A_cell')
## END


d[["subcluster"]]<-Idents(d)

meta<-d@meta.data

meta$subcluster<-as.character(meta$subcluster)

## annotate cluster 0
clusters<-c(1)
for (i in clusters){
  meta_sub<-readRDS(paste0("Cluster_",i,"/metadata_C.rds"))
  subcluster_type<-unique(meta_sub$subcluster)
  for (st in subcluster_type){
    cells<-rownames(meta_sub)[which(meta_sub$subcluster==st)]
	meta$subcluster[which(rownames(meta) %in% cells)]<-st
  }
}

d[["subcluster"]]<-meta$subcluster
Idents(d)<-d[["subcluster"]]

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



levels(d) <- c('A_cell','D_cell','G_cell','L_cell','M_cell','EC_cell','EC_cell_like')

metadata<-d@meta.data
enriched_idents(metadata,idents_tested="Histological",idents_based="subcluster",tmpdir="enriched_idents_Histological_subcluster",angle=90)
enriched_idents(metadata,idents_tested="subcluster",idents_based="Histological",tmpdir="enriched_idents_subcluster_Histological",angle=90,width=10)


features<-c('CHGA','CHGB','GHRL','SST','GAST','GCG','MLN','TPH1','DDC','HDC')


features<-unique(c(features))

g<-DotPlot(d, 
      features = features, 
      cols = c('lightgrey', 'red')) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
ggsave("DotPlot_markers_annotated.png",width=10,height=6,plot=g)





dir.create("DEG_subclustered")
setwd("DEG_subclustered")

Idents(d)<-d[["Histological"]]
sc_all<-unique(d[["subcluster"]]$subcluster)

for (sc in sc_all){
  print(paste0("Cell type: ",sc))
  dd<-subset(d,subcluster==sc)
  
  dir.create(sc)
  
## using logCPM normalization data; wilcox method
  celltype.markers<-FindAllMarkers(dd,assay="RNA",slot="data",max.cells.per.ident=Inf,only.pos=FALSE)

  celltype.markers<-celltype.markers[,c("cluster","avg_log2FC","pct.1","pct.2","p_val","p_val_adj","gene")]
  celltype.markers<-celltype.markers[which(celltype.markers$p_val_adj<0.05),]

  celltype.markers.neg<-subset(celltype.markers,avg_log2FC<0)
  celltype.markers<-subset(celltype.markers,avg_log2FC>0)

  celltype.markers.neg<-celltype.markers.neg[order(celltype.markers.neg$cluster,celltype.markers.neg$avg_log2FC),]
  write.table(celltype.markers.neg,paste0(sc,"/celltype_NEG_markers_Wilcox_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)

  celltype.markers<-celltype.markers[order(celltype.markers$cluster,-celltype.markers$avg_log2FC),]
  
  write.table(celltype.markers,paste0(sc,"/celltype_markers_Wilcox_all.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  celltype.markers.top20<-NULL
  for (j in unique(celltype.markers$cluster)){
    k<-which(celltype.markers$cluster==j)
    if (length(k)>20){
      k<-k[c(1:20)]
    }
    celltype.markers.top20<-rbind(celltype.markers.top20,celltype.markers[k,])
  }
  write.table(celltype.markers.top20,paste0(sc,"/celltype_markers_Wilcox_top20.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  
  rm(dd)
  gc()
}

