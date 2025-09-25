
norm01<-function(x,mini=NULL,maxi=NULL){
  if (min(x)<max(x)){
    if (is.null(mini)){
      mini<-min(x)
    }
    if (is.null(maxi)){
      maxi<-max(x)
    }
    return((x-mini)/(maxi-mini))
  }else{
    return(rep(0.5,length(x)))
  }
}



DEgeneOrderSimilarity<-function(query,reference){
  ## mapping query cluster expression markers to a reference (known) cell-type marker list and calculate a similarity score
  ## Jaccard Coefficient Like method
  score_list<-numeric(length(reference))
  for (i in 1:length(reference)){
    ref_sub<-reference[c(1:i)]
    score<-length(intersect(query,ref_sub))/length(query)*(1-(i/length(reference)))
    score_list[i]<-score
  }
  return(score_list)
}



Heatmap.Downsample<-function(seurat,max.cell.per.cluster=1000,min.cell.per.cluster=50,remove.unassigned=TRUE,random.seed=1,DSbyRatio=TRUE){
  if (!is.na(random.seed)){
    set.seed(random.seed)
  }
  require(Seurat)
  cellnames<-NULL
  allnames<-CellsByIdentities(object = seurat)
  #cluster.names<-names(allnames)
  cluster.cells.numbers<-sort(table(Idents(seurat)),decreasing=TRUE)
  max.cells<-0
  for (i in 1:length(allnames)){
    if (length(allnames[[i]])>max.cells){
      max.cells<-length(allnames[[i]])
    }
  }
  if (max.cells>max.cell.per.cluster){
    downsample.ratio<-max.cell.per.cluster/max.cells
    for (i in names(cluster.cells.numbers)){
      if (remove.unassigned & i=="unassigned"){
        next
      }
      cluster.cells<-allnames[[i]]
      num<-trunc(downsample.ratio*length(cluster.cells))
      if (num<min.cell.per.cluster){
        num<-min.cell.per.cluster
      }
      if (!DSbyRatio & num<max.cell.per.cluster){ ## by default down-sample according to cell proportion in the cluster, by can also be set to down-sample to max.cell.per.cluster whenever possible
        num<-max.cell.per.cluster
      }
      if (num>length(cluster.cells)){
        num<-length(cluster.cells)
      }
      cellnames<-c(cellnames,cluster.cells[sample(c(1:length(cluster.cells)),num,replace=FALSE)])
      print(paste0("Cell type ",i," No.: ",num))
    }
  }else{
    print("Down-sample not applied.")
    cellnames<-Cells(seurat)
  }
  return(cellnames)
}



Heatmap.Downsample.Simple<-function(seurat,max.cell.per.cluster=1000,min.cell.per.cluster=50,remove.unassigned=TRUE,random.seed=1,DSbyRatio=TRUE){
  ## simple version of Heatmap.Downsample
  ## input a data.frame with two cols: cells, idents
  
  if (!is.na(random.seed)){
    set.seed(random.seed)
  }
  
  cellnames<-NULL
  # allnames<-seurat$cells
  cluster.names<-unique(seurat$idents)
  cluster.cells.numbers<-sort(table(seurat$idents),decreasing=TRUE)
  # max.cells<-0
  # for (i in 1:length(allnames)){
    # if (length(allnames[[i]])>max.cells){
      # max.cells<-length(allnames[[i]])
    # }
  # }
  max.cells<-cluster.cells.numbers[1]
  # print(max.cells)
  if (max.cells>max.cell.per.cluster){
    downsample.ratio<-max.cell.per.cluster/max.cells
    for (i in names(cluster.cells.numbers)){
      if (remove.unassigned & i=="unassigned"){
        next
      }
      # cluster.cells<-allnames[[i]]
      cluster.cells<-seurat$cells[which(seurat$idents==i)]
      num<-trunc(downsample.ratio*length(cluster.cells))
      if (num<min.cell.per.cluster){
        num<-min.cell.per.cluster
      }
      if (!DSbyRatio & num<max.cell.per.cluster){ ## by default down-sample according to cell proportion in the cluster, by can also be set to down-sample to max.cell.per.cluster whenever possible
        num<-max.cell.per.cluster
      }
      if (num>length(cluster.cells)){
        num<-length(cluster.cells)
      }
      # print(length(cluster.cells))
      # print(num)
      cellnames<-c(cellnames,cluster.cells[sample(c(1:length(cluster.cells)),num,replace=FALSE)])
      # print(paste0("Cell type ",i," No.: ",num))
    }
  }else{
    # print("Down-sample not applied.")
    cellnames<-seurat$cells
  }
  return(cellnames)
}



mergeSameGene_avg<-function(exp_matrix,genes){
  ## input an expression matrix with samples as cols and Genes as rows; values are directly averaged for rows with the same gene name
  dups<-union(which(duplicated(genes,fromLast=TRUE)),which(duplicated(genes)))
  uni<-setdiff(c(1:nrow(exp_matrix)),dups)
  uni_matrix<-exp_matrix[uni,]
  rownames(uni_matrix)<-genes[uni]
  dup_genes<-unique(genes[dups])
  dup_matrix<-matrix(ncol=ncol(exp_matrix),nrow=length(dup_genes))
  rownames(dup_matrix)<-dup_genes
  colnames(dup_matrix)<-colnames(exp_matrix)
  for (gene in dup_genes){
    tmp<-exp_matrix[which(genes==gene),]
    dup_matrix[gene,]<-apply(tmp,2,mean)
  }
  exp_matrix_new<-rbind(uni_matrix,dup_matrix)
  return(exp_matrix_new)

}




mergeSameGene_sum<-function(exp_matrix,genes){
  ## input an expression matrix with samples as cols and Genes as rows; similar with function "mergeSameGene_avg", 
  ## BUT values are directly summed instead of averaged for rows with the same gene name
  dups<-union(which(duplicated(genes,fromLast=TRUE)),which(duplicated(genes)))
  uni<-setdiff(c(1:nrow(exp_matrix)),dups)
  uni_matrix<-exp_matrix[uni,]
  rownames(uni_matrix)<-genes[uni]
  dup_genes<-unique(genes[dups])
  dup_matrix<-matrix(ncol=ncol(exp_matrix),nrow=length(dup_genes))
  rownames(dup_matrix)<-dup_genes
  colnames(dup_matrix)<-colnames(exp_matrix)
  for (gene in dup_genes){
    tmp<-exp_matrix[which(genes==gene),]
    dup_matrix[gene,]<-apply(tmp,2,sum)
  }
  exp_matrix_new<-rbind(uni_matrix,dup_matrix)
  return(exp_matrix_new)
}



ensemblID2geneSymbol<-function(query,reference,species="human"){
  ## query: a character vector containing gene ensembl IDs
  ## reference: a data.frame with cols "ensemblID","geneSymbol"
  if (species=="human"){
    query<-substr(query,1,15)  ## standard human gene ensembl IDs length
  }else{
    if (species=="mouse"){
      query<-substr(query,1,18)  ## standard mouse gene ensembl IDs length
    }
  }
  ref<-reference$geneSymbol
  names(ref)<-reference$ensemblID
  res<-ref[query]
  names(res)<-NULL
  return(res)
}

makeEnsemblID2GeneSymbolFromGTF<-function(gtf,format="gtf",filename="ensemblID2geneSymbol.txt"){
  library(rtracklayer)
  tab<-import(gtf,format=format)
  tab<-as.data.frame(tab)
  tab<-subset(tab,type=="gene")
  tab<-tab[,c("gene_id","gene_name")]
  tab<-unique(tab)
  colnames(tab)<-c("ensemblID","geneSymbol")
  if (!is.na(filename)){
    write.table(tab,filename,sep="\t",row.names=FALSE,quote=FALSE)
  }
  return(tab)
}


CPM<-function(dx){
  if (is.matrix(dx)){
    mat<-TRUE
    dz<-data.frame(dx)
    geneSymbol<-row.names(dx)
    colnames(dz)<-colnames(dx)
    dz<-cbind(geneSymbol,dz)
    row.names(dz)<-NULL
    dx<-dz
    rm(dz)
    gc()
  }else{
    mat<-FALSE
  }

  for (i in 2:ncol(dx)){
    dx[,i]<-dx[,i]/sum(dx[,i])*1e6
  }
  if (mat){
    dz<-as.matrix(dx[,2:ncol(dx)])
    colnames(dz)<-colnames(dx)[2:ncol(dx)]
    row.names(dz)<-dx$geneSymbol
    return(dz)
  }else{
    return(dx)
  }
}


TPM<-function(dx,geneLength){
  if (is.matrix(dx)){
    mat<-TRUE
    dz<-data.frame(dx)
    geneSymbol<-row.names(dx)
    colnames(dz)<-colnames(dx)
    dz<-cbind(geneSymbol,dz)
    row.names(dz)<-NULL
    dx<-dz
    rm(dz)
    gc()
  }else{
    mat<-FALSE
  }
  if (!identical(dx$geneSymbol,geneLength$geneSymbol)){
    print("WARNING: gene name order or number not correct! Try to correct it now...")
    q<-union(setdiff(dx$geneSymbol,geneLength$geneSymbol),setdiff(geneLength$geneSymbol,dx$geneSymbol))
    if (length(q)>0){
      print("Some genes not found in both files are excluded first...")
      z<-integer(nrow(geneLength))
      names(z)<-geneLength$geneSymbol
      z[intersect(q,geneLength$geneSymbol)]<-1
      geneLength<-geneLength[which(z==0),]
      
      z<-integer(nrow(dx))
      names(z)<-dx$geneSymbol
      z[intersect(q,dx$geneSymbol)]<-1
      dx<-dx[which(z==0),]
    }
    geneLength$geneSymbol<-factor(geneLength$geneSymbol,levels=dx$geneSymbol,ordered=T)
    geneLength<-geneLength[order(geneLength$geneSymbol),]
  }
  for (i in 2:ncol(dx)){
    dx[,i]<-dx[,i]/geneLength$geneLength
    dx[,i]<-dx[,i]/sum(dx[,i])*1e6
  }
  if (mat){
    dz<-as.matrix(dx[,2:ncol(dx)])
    colnames(dz)<-colnames(dx)[2:ncol(dx)]
    row.names(dz)<-dx$geneSymbol
    return(dz)
  }else{
    return(dx)
  }
}


TPM_single<-function(x,geneLength){
  ## NO gene name or order correspondence check
  x<-x/geneLength$geneLength
  x<-x/sum(x)*1e6
  return(x)
}



enriched_idents<-function(metadata,idents_tested,idents_based,idents_tested_order=NULL,idents_based_order=NULL,
                          do_plot=TRUE,tmpdir="enriched_idents_tmp",plot_device="png",width=7,height=7,unit="in",
                          test_adjust_method="bonferroni",test_method="chisq",no_marks=FALSE,
                          high_color='#f16a27',low_color='#ffffff',midpoint=NULL,mid_color=NULL,balanced.scale.bar=FALSE,logROE=FALSE,angle=0){
  # metadata: meta-data of the seurat object
  # idents_tested: the feature to be tested for enrichment
  # idents_based: the feature that split cells into different groups
  
  ## the rationale for using Fisher's exact test for cell type enrichment test:
  ## ref: DOI: 10.1038/s41556-020-00613-6; Extended Fig 3d (right) "Cluster 12 showed significant enrichment (Fisher’s exact test, P = 6.7 × 10−41 odds ratio (OR) = 5.7) for cells in G2/M (Extended Data Fig. 3d)", "contingency table showing number of cells in G2/M for CellCycle subtype compared to all the other cells (right), One-sided P value calculated with Fisher’s exact test"
  ## ref: https://www.ccjm.org/content/84/9_suppl_2/e20 ; https://www.datascienceblog.net/post/statistical_test/contingency_table_tests/
  ## ref: DOI: 10.1038/s41422-022-00627-9 ; in this article the R(O/E) of cell type distribution was calculated by chi-squared test
  ## ref: DOI: 10.1016/j.cell.2020.03.048; also use chi-squared test
  
  require(data.table)
  require(ggplot2)
  require(ggsci)
  if (!dir.exists(tmpdir)){
    dir.create(tmpdir)
  }
  y<-c(idents_tested,idents_based)
  metadata<-as.data.table(metadata)
  metadata<-metadata[,..y]
  colnames(metadata)<-c("idents_tested","idents_based")
  print(table(metadata$idents_based))
  print(table(metadata$idents_tested))
  metadata$idents_tested<-as.character(metadata$idents_tested)
  metadata$idents_based<-as.character(metadata$idents_based)
  metadata<-subset(metadata,!is.na(idents_tested))
  metadata<-subset(metadata,!is.na(idents_based))
  summary_stats<-data.table(idents_tested=unique(metadata$idents_tested))
  setkey(summary_stats,idents_tested)
  z<-unique(metadata$idents_based)
  for (i in 1:length(z)){
    s<-as.data.table(table(metadata[idents_based==z[i],idents_tested]))
    colnames(s)<-c("idents_tested",z[i])
    s[is.na(idents_tested),2]
    setkey(s,idents_tested)
    summary_stats<-merge(summary_stats,s,all.x=TRUE)
  }
  
  setnafill(summary_stats, fill = 0,cols=c(2:ncol(summary_stats)))
  ## apply fisher-exact test (or chi-squared test)
  res<-NULL
  summary_stats_fisher<-as.matrix(summary_stats,rownames="idents_tested")
  cs<-colSums(summary_stats_fisher)
  names(cs)<-colnames(summary_stats_fisher)
  rs<-rowSums(summary_stats_fisher)
  names(rs)<-rownames(summary_stats_fisher)
  su<-sum(summary_stats_fisher)
  
  g2<-NULL
  
  if (test_method=="fisher"){
    for (i in rownames(summary_stats_fisher)){
      for (j in colnames(summary_stats_fisher)){
        res_matrix<-matrix(ncol=2,nrow=2)
        rownames(res_matrix)<-c(i,"idents_based_other")
        colnames(res_matrix)<-c(j,"idents_tested_other")
        res_matrix[i,j]<-summary_stats_fisher[i,j]
        res_matrix[i,"idents_tested_other"]<-rs[i]-summary_stats_fisher[i,j]
        res_matrix["idents_based_other",j]<-cs[j]-summary_stats_fisher[i,j]
        res_matrix["idents_based_other","idents_tested_other"]<-su-cs[j]-rs[i]+summary_stats_fisher[i,j]
        OR<-res_matrix[1,1]/res_matrix[1,2]/res_matrix[2,1]*res_matrix[2,2]
        test<-fisher.test(res_matrix,alternative="greater")
        res_sub<-data.frame(idents_based=i,idents_tested=j,ingroup=res_matrix[1,1],outgroup=res_matrix[2,1],ingroup_others=res_matrix[1,2],outgroup_others=res_matrix[2,2],OR=OR,pval=test$p.value)

        res<-rbind(res,res_sub)
      }
    }
    res$padj<-p.adjust(res$pval,method=test_adjust_method)
    res$sig_mark<-""

    sig_mark<-c("*","**","***")
    sig_base<-c(5e-2,1e-2,1e-3)
    names(sig_mark)<-sig_base
    for (i in sig_base){
      res$sig_mark[which(res$padj<i)]<-sig_mark[as.character(i)]
    }
  }else{
    ## Chi-squared test
    test<-chisq.test(summary_stats_fisher)
    
    ## idents_tested as Rows and idents_based as Cols
    
    print(test)
    ROE<-test$observed/test$expected

    melt_obs<-reshape2::melt(test$observed)
    melt_exp<-reshape2::melt(test$expected)
    melt_roe<-reshape2::melt(ROE)
    
    melt_merge<-melt_obs
    colnames(melt_merge)<-c('idents_tested','idents_based','observed')
    melt_merge$expected<-melt_exp$value
    melt_merge$ROE<-melt_roe$value
    
    res<-melt_merge
    
    res$enrichment_mark<-"-"
    
    z<-which(res$ROE>0 & res$ROE<0.2)
    res$enrichment_mark[z]<-"+/-"
    
    z<-which(res$ROE>=0.2 & res$ROE<=0.8)
    res$enrichment_mark[z]<-"+"
    
    z<-which(res$ROE>0.8 & res$ROE<=1)
    res$enrichment_mark[z]<-"++"
    
    z<-which(res$ROE>1)
    res$enrichment_mark[z]<-"+++"
    
    if (!is.null(idents_tested_order)){
      res$idents_tested<-factor(res$idents_tested,levels=idents_tested_order,ordered=TRUE)
    }
    if (!is.null(idents_based_order)){
      res$idents_based<-factor(res$idents_based,levels=idents_based_order,ordered=TRUE)
    }
    
    res$logROE<-log2(res$ROE+1e-2)
    
    if (! logROE){
      g2<-ggplot(res)+
        geom_tile(aes(x=idents_based,y=idents_tested,fill=ROE))
        values<-res$ROE
        legend_name<-'Ro/e'
    }else{
      g2<-ggplot(res)+
        geom_tile(aes(x=idents_based,y=idents_tested,fill=logROE))
        values<-res$logROE
        legend_name<-'log2 Ro/e'
    }
    
    if (!no_marks){
      g2<-g2+geom_text(aes(x=idents_based,y=idents_tested,label=enrichment_mark))
    }
    
    if (is.null(mid_color)){
      g2<-g2+scale_fill_gradient(high=high_color,low=low_color,name=legend_name)
    }else{
      g2<-g2+scale_fill_gradient2(high=high_color,low=low_color,mid=mid_color,midpoint=midpoint,name=legend_name)
    }
    
    g2<-g2+theme_classic()+xlab(idents_based)+ylab(idents_tested)
  }
  
  
  ## END
  
  summary_stats_prop<-as.matrix(summary_stats,rownames="idents_tested")
  summary_stats_prop<-apply(summary_stats_prop,2,function(x){return(x/sum(x))})
  summary_stats_prop<-summary_stats_prop*100
  
  
  summary_stats_prop<-reshape2::melt(summary_stats_prop)
  colnames(summary_stats_prop)<-c("idents_tested","idents_based","Counts")
  summary_stats_prop$idents_based<-factor(summary_stats_prop$idents_based)
  summary_stats_prop$idents_tested<-factor(summary_stats_prop$idents_tested)

  g<-ggplot(summary_stats_prop)+geom_bar(aes(x=idents_based,weight=Counts,fill=idents_tested),position="stack")+
     theme(axis.text=element_text(size=13),axis.text.x=element_text(size=13,angle=angle,vjust=0.5),axis.title=element_text(size=16),legend.text=element_text(size=13),legend.title=element_text(size=16))+
     ylab("Proportion (%)")+xlab(idents_based)+
     theme_classic()
     
  if (length(levels(summary_stats_prop$idents_tested))>26){
    g<-g+scale_fill_discrete(name=idents_tested)
  }else{
    g<-g+scale_fill_manual(values=pal_ucscgb("default")(26),name=idents_tested)
  }
     

  ggsave(filename=paste0(tmpdir,"/stack_barplot.",plot_device),plot=g,device=plot_device,width=width,height=height,unit=unit)
  if (test_method=="fisher"){
    write.table(res,paste0(tmpdir,"/fisher_test.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  }else{
    write.table(res,paste0(tmpdir,"/chisq.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  }
  
  return(list(res=res,g=g,g2=g2))
}



plot_heatmap_avg_exps<-function(exps,class_label,select_genes=NULL,cluster_order=NULL,gene_order=NULL,scaling=TRUE,centering=TRUE,ignore_class_sample_size=TRUE,
                                high_col='#f08c10',low_col='#131753',mid_col='#ffffff',disable_auto_theme=FALSE){
  require(ggplot2)
  if (!is.null(select_genes)){
    exps<-exps[select_genes,]
  }
  
  if (centering & !ignore_class_sample_size){
    res<-generate_centroid_by_class(exps=exps,class_label=class_label,scaling=scaling,centering=centering)
  }else{
    res<-generate_centroid_by_class(exps=exps,class_label=class_label,scaling=FALSE,centering=FALSE)
  }
  
  avg_exps<-res$centroids
  
  ## centering and scaling after averaging (ignoring the differnet sample size of each class)
  if (centering & ignore_class_sample_size){
    avg_exps<-t(scale(t(avg_exps),center=centering,scale=scaling))
  }
  
  exps_melt<-reshape2::melt(avg_exps)
  
  colnames(exps_melt)<-c('gene','cluster','exps')
  
  if (!is.null(cluster_order)){
    exps_melt$cluster<-factor(exps_melt$cluster,levels=cluster_order,ordered=TRUE)
  }
  
  if (!is.null(gene_order)){  ## pre-define row (gene) order
    exps_melt$gene<-factor(exps_melt$gene,levels=gene_order,ordered=TRUE)
  }
  
  md<-(max(exps_melt$exps)-min(exps_melt$exps))/2+min(exps_melt$exps)
  
  g<-ggplot(exps_melt)+
     geom_tile(aes(x=cluster,y=gene,fill=exps))+
     scale_fill_gradient2(high=high_col,low=low_col,mid=mid_col,midpoint=md)
  if (!disable_auto_theme){
    g<-g+scale_x_discrete(position = "top")+
       theme_classic()+
       theme(axis.text.x.top=element_text(angle=90,vjust=0.5))
  }
  return(list(avg_exps=avg_exps,g=g))
}



DE_malignant_method<-function(seurat.obj,bulk_normal_markers=NULL,bulk_tumor_markers=NULL,num_markers=100,normal_cells=NULL,padj=0.01,log2FC=NA,tmpdir="DE_malignant_tmp",
                              slot="data",assay="RNA",force_normal_orig_as_nonMag=TRUE,markers_tolerate=0,cell_tolerate=0,malignant_marker_constraint=NULL,non_malignant_marker_constraint=NULL,scaled01=FALSE){
  ## ref: DE method from PMID:32532891
  ## modified: use normal tissue cells as reference so bulk markers for the initial run not necessary
  ## modified: add constraints for marker exploration to those already elucinated in bulk study (as custom provided)
  
  norm01<-function(x){
    if (min(x)<max(x)){
      return((x-min(x))/(max(x)-min(x)))
    }else{
      return(rep(0.5,length(x)))
    }
  }
  
  malignant_cell<-NULL
  
  if (!is.null(malignant_marker_constraint) & !is.null(non_malignant_marker_constraint)){
    marker_constraint<-TRUE
  }else{
    marker_constraint<-FALSE
  }
  
  require(Seurat)
  require(ggplot2)
  print("Start DE-malignant-method analysis...")
  if (!dir.exists(tmpdir)){
    dir.create(tmpdir)
  }
  DefaultAssay(seurat.obj)<-assay  ## switch back to un-corrected data

  if (!is.null(bulk_normal_markers) & !is.null(bulk_tumor_markers)){
    if (is.na(num_markers)){
      num_markers<-length(bulk_tumor_markers)
    }else{
      if (num_markers<length(bulk_tumor_markers)){
        bulk_tumor_markers<-bulk_tumor_markers[1:num_markers]
      }
      if (num_markers<length(bulk_normal_markers)){
        bulk_normal_markers<-bulk_normal_markers[1:num_markers]
      }
    }
    ## initial run: use TCGA bulk DE markers
    normal_markers <- bulk_normal_markers
    tumor_markers<- bulk_tumor_markers
  }else{
    ## no initial bulk marker provided
    print("Using cells from normal tissues for initial marker set detection")
    normal_markers<-NULL
    tumor_markers<-NULL
  }
  
  stable<-FALSE
  run<-0
  while(!stable){
    run<-run+1
    print(paste0("Run: ",run))
    
    if (!dir.exists(paste0(tmpdir,"/run_",run))){
      dir.create(paste0(tmpdir,"/run_",run))
    }
    
    if (!is.null(normal_markers) & !is.null(tumor_markers)){
      
      seurat.obj <- AddModuleScore(
      object = seurat.obj,
        features = list(normal_markers),
        name = 'non_malignant_score'
      )

      seurat.obj <- AddModuleScore(
        object = seurat.obj,
        features = list(tumor_markers),
        name = 'malignant_score'
      )
      png(paste0(tmpdir,"/run_",run,"/UMAP_Non_Mag_Score.png"),width=800,height=800)
      print(FeaturePlot(object = seurat.obj, features = 'non_malignant_score1',reduction="umap"))
      dev.off()
      png(paste0(tmpdir,"/run_",run,"/UMAP_Mag_Score.png"),width=800,height=800)
      print(FeaturePlot(object = seurat.obj, features = 'malignant_score1',reduction="umap"))
      dev.off()
      png(paste0(tmpdir,"/run_",run,"/tSNE_Non_Mag_Score.png"),width=800,height=800)
      print(FeaturePlot(object = seurat.obj, features = 'non_malignant_score1',reduction="tsne"))
      dev.off()
      png(paste0(tmpdir,"/run_",run,"/tSNE_Mag_Score.png"),width=800,height=800)
      print(FeaturePlot(object = seurat.obj, features = 'malignant_score1',reduction="tsne"))
      dev.off()
      res<-data.frame(non_malignant_score1=seurat.obj[["non_malignant_score1"]],malignant_score1=seurat.obj[["malignant_score1"]])
      if (scaled01){
        res$non_malignant_score1<-norm01(res$non_malignant_score1)
        res$malignant_score1<-norm01(res$malignant_score1)
      }
      rownames(res)<-colnames(seurat.obj)
      km<-kmeans(res,centers=2,nstart=20)
      means<-as.data.frame(km$centers)
      print(means)
      if (((means[2,1]-means[1,1])*(means[2,2]-means[1,2]))>0){
        warning("Clusters may NOT be distinguished.")
        print("WARNING: Clusters may NOT be distinguished.")
      }
      ## classified the cluster with higher malignant_score as malignant
      mal<-which(means$malignant_score1>=max(means$malignant_score1))
      res$cluster_num<-km$cluster
      res$cluster<-"non_malignant"
      res$cluster[which(res$cluster_num==mal)]<-"malignant"
      if (force_normal_orig_as_nonMag & !is.null(normal_cells)){
        res$cluster[which(rownames(res) %in% normal_cells)]<-"non_malignant"
      }
      res$cluster<-factor(res$cluster)
      res$cluster<-relevel(res$cluster,ref="non_malignant")

      g<-ggplot(res)+geom_point(aes(x=non_malignant_score1,y=malignant_score1,color=cluster),shape=16)+
      geom_point(data=means,mapping=aes(x=non_malignant_score1,y=malignant_score1),color="black",shape=9,size=5)+
      theme(axis.text=element_text(size=10),axis.title=element_text(size=15),legend.text=element_text(size=12),legend.title=element_text(size=17))
      
      png(paste0(tmpdir,"/run_",run,"/Score_pointPlt.png"),width=800,height=800)
      print(g)
      dev.off()
    }else{
      ## use normal cells as reference
      res<-data.frame(cluster=rep("malignant",ncol(seurat.obj)))
      rownames(res)<-colnames(seurat.obj)
      res$cluster[which(rownames(res) %in% normal_cells)]<-"non_malignant"
      res$cluster<-factor(res$cluster)
      res$cluster<-relevel(res$cluster,ref="non_malignant")
    }
    Idents(seurat.obj)<-res$cluster
    
    malignant_cell_now<-which(res$cluster=="malignant")
    
    png(paste0(tmpdir,"/run_",run,"/tSNE_current_definition.png"),width=800,height=800)
    print(DimPlot(object = seurat.obj,reduction="tsne"))
    dev.off()
    
    png(paste0(tmpdir,"/run_",run,"/UMAP_current_definition.png"),width=800,height=800)
    print(DimPlot(object = seurat.obj,reduction="umap"))
    dev.off()
    
    markers_all<-FindMarkers(seurat.obj,ident.1="malignant",only.pos=FALSE)  ## use non_malignant as reference; avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
    markers_all<-subset(markers_all,p_val_adj<padj)
    markers_all<-markers_all[order(abs(markers_all$avg_log2FC),decreasing=TRUE),]
    if (!is.na(log2FC)){
      markers_all<-subset(markers_all,abs(avg_log2FC)>=log2FC)
    }
    fit.res.tumorHE<-markers_all[which(markers_all$avg_log2FC>0),]
    fit.res.normalHE<-markers_all[which(markers_all$avg_log2FC<0),]
    if (!is.na(num_markers)){
      normal_markers_new<-rownames(fit.res.normalHE)[1:num_markers]
      tumor_markers_new<-rownames(fit.res.tumorHE)[1:num_markers]
    }else{
      normal_markers_new<-rownames(fit.res.normalHE)
      tumor_markers_new<-rownames(fit.res.tumorHE)
    }
    
    normal_markers_new<-normal_markers_new[which(!is.na(normal_markers_new))]  # remove NA
    tumor_markers_new<-tumor_markers_new[which(!is.na(tumor_markers_new))]  # remove NA
    
    if (marker_constraint){
      normal_markers_new<-intersect(normal_markers_new,non_malignant_marker_constraint)
      tumor_markers_new<-intersect(tumor_markers_new,malignant_marker_constraint)
    }
    
    diff_normal<-c(setdiff(normal_markers_new,normal_markers),setdiff(normal_markers,normal_markers_new))
    diff_tumor<-c(setdiff(tumor_markers_new,tumor_markers),setdiff(tumor_markers,tumor_markers_new))
    
    if (((length(diff_normal)+length(diff_tumor))<=markers_tolerate) & ((length(setdiff(malignant_cell_now,malignant_cell))+length(setdiff(malignant_cell,malignant_cell_now)))<=cell_tolerate)){
      print("Result stabled.")
      stable<-TRUE
      if (((means[2,1]-means[1,1])*(means[2,2]-means[1,2]))>0){
        print("WARNING: Results may not be desired. The cluster with higher malignant score also has higher non-malignant score!")
      }
    }else{
      print("New malignant markers:")
      print(tumor_markers_new)
      #print(paste0("New malignant markers No.: ",length(tumor_markers_new)))
      print(paste0("Different malignant markers No.: ",length(diff_tumor)))
      print("New non-malignant markers:")
      print(normal_markers_new)
      #print(paste0("New non-malignant markers No.: ",length(normal_markers_new)))
      print(paste0("Different non-malignant markers No.: ",length(diff_normal)))
      print(paste0("Malignant cell No.: ",length(malignant_cell_now)))
    }
    
    malignant_cell<-malignant_cell_now
    normal_markers<-normal_markers_new
    tumor_markers<-tumor_markers_new
    cat("",sep="\n")
    if (length(normal_markers)!=length(tumor_markers)){
      imp<-rep("",abs(length(normal_markers)-length(tumor_markers)))
      if (length(normal_markers)<length(tumor_markers)){
        current_markers<-data.frame(normal_markers=c(normal_markers,imp),tumor_markers=tumor_markers)
      }else{
        current_markers<-data.frame(normal_markers=normal_markers,tumor_markers=c(tumor_markers,imp))
      }
    }else{
      current_markers<-data.frame(normal_markers=normal_markers,tumor_markers=tumor_markers)
    }
    
    write.table(current_markers,paste0(tmpdir,"/run_",run,"/current_markers.txt"),sep="\t",row.names=FALSE,quote=FALSE)
  }
  if (!is.null(bulk_normal_markers) & !is.null(bulk_tumor_markers)){
    print("Intersection between obtained normal markers and original Bulk normal markers:")
    print(intersect(normal_markers,bulk_normal_markers))
    print(paste0("Intersection Ratio: ",round(length(intersect(normal_markers,bulk_normal_markers))/length(normal_markers)*100,digits=2),"%"))
    print("Intersection between obtained tumor markers and original Bulk tumor markers:")
    print(intersect(tumor_markers,bulk_tumor_markers))
    print(paste0("Intersection Ratio: ",round(length(intersect(tumor_markers,bulk_tumor_markers))/length(tumor_markers)*100,digits=2),"%"))
    
    gc()
  }
  return(list(cellType.res=res,normal_markers=normal_markers,tumor_markers=tumor_markers))
}



scCancer_malignant_score<-function(CNVObj,ncores=4,method="pearson",n_neighbors=100,gene_chr_annotation="/data3/liangjq/scRNA-seq/gen_pos.txt"){
  if (ncores>1){
    require(parallel)
  }
  require(WGCNA)
  infercnv_exp<-CNVObj@expr.data  ## rows as genes and cols as cells
  rm(CNVObj)
  gc()
  
  print("Calculating Similarity Matrix...")
  Similarities<-WGCNA::cor(infercnv_exp,nThreads=ncores,method=method)
  diag(Similarities)<-0  ## mask self-correlation
  
  gc()
  
  boundary_mask_scaled<-function(x){
    s<-Similarities[,x]
    bound<-sort(s,decreasing=TRUE)[(n_neighbors+1)]  ## similaries < similaries[n_neighbors] are masked
    s[which(s<=bound)]<-0
    s<-s/sum(s)  ## re-scale to sum to 1
    return(s)
  }
  print("Similarity Matrix mask and scaling...")
  if (ncores>1){
    Similarities_res<-mclapply(c(1:ncol(Similarities)),boundary_mask_scaled,mc.preschedule=TRUE,mc.cores=ncores)
  }else{
    Similarities_res<-lapply(c(1:ncol(Similarities)),boundary_mask_scaled)
  }
  Similarities_res<-unlist(Similarities_res)
  Similarities<-matrix(Similarities_res,ncol=ncol(Similarities))
  
  rm(Similarities_res)
  gc()
  
  diag(Similarities)<-1  ## restore self-correlation as 1
  
  infercnv_exp_new<-matrix(0,ncol=ncol(infercnv_exp),nrow=nrow(infercnv_exp),
                           dimnames=list(rownames(infercnv_exp),colnames(infercnv_exp)))
  
  print("Calculating averaged CNV value for cells")
  for (i in 1:ncol(infercnv_exp_new)){
    infercnv_exp_new[,i]<-rowSums(infercnv_exp*Similarities[i,])
  }
  
  gene_chr<-read.table(gene_chr_annotation,header=FALSE,sep="\t")
  gene_chr<-gene_chr[,c(1,2)]
  colnames(gene_chr)<-c("gene","chr")
  gene_chr<-gene_chr[which(gene_chr$gene %in% rownames(infercnv_exp_new)),]
  gene_chr$gene<-factor(gene_chr$gene,levels=rownames(infercnv_exp_new),ordered=TRUE)
  gene_chr<-gene_chr[order(gene_chr$gene),]
  unique_chr<-unique(as.character(gene_chr$chr))
  
  CNV_per_chr<-matrix(0,nrow=length(unique_chr),ncol=colnames(infercnv_exp_new),
                      dimnames=list(unique_chr,colnames(infercnv_exp_new)))
  
  for (chr in unique_chr){
    genes_in_chr<-which(gene_chr$chr==chr)
    CNV_per_chr[chr,]<-apply(infercnv_exp_new[genes_in_chr,],2,var)
  }
  
  print("Calculating final malignant scores...")
  score<-colMeans(CNV_per_chr)
  names(score)<-colnames(infercnv_exp_new)
  
  return(list(malignant_score=score,averaged_CNV_values=infercnv_exp_new,averaged_CNV_perCHR=CNV_per_chr))
}


malignant_score_custom<-function(cna_matrix,cell_idents,x.range=c(-Inf,Inf),score.range=c(0,NA),plot_width=8,plot_height=6,angle=45,filename="malignant_score",violin_adjust=1){
  require(ggplot2)
  #require(ggsci)

  if (is.data.frame(cell_idents)){
    cell_idents<-unlist(cell_idents)
  }
  
  if (is.finite(x.range[1]) | is.finite(x.range[2])){
  print(paste0("Clipped values outside x.range (",x.range[1],",",x.range[2],")"))
  cna_matrix[cna_matrix<x.range[1]]<-x.range[1]
  cna_matrix[cna_matrix>x.range[2]]<-x.range[2]
  }
  
  med<-median(cna_matrix)
  cna_matrix<-cna_matrix-med  ## substract median first
  
  score<-apply(cna_matrix,2,var)  ## calculate per-cell CNV variance
  
  z<-unique(cell_idents)
  score_med<-data.frame(cell_idents=z,score=numeric(length(z)))
  for (i in 1:length(z)){
    score_med$score[i]<-median(score[which(cell_idents==z[i])])
  }
  print(score_med)
  
  ord<-score_med$cell_idents[order(score_med$score,decreasing=TRUE)]
  
  dataGG<-data.frame(cell=colnames(cna_matrix),score=score,identity=factor(cell_idents,levels=ord,ordered=TRUE))
  g1<-ggplot(dataGG)+
     geom_density(aes(x=score,fill=identity),alpha=0.5)+
     geom_vline(aes(xintercept=score,color=cell_idents),data=score_med,size=0.8,linetype="longdash")+
     theme(axis.text=element_text(size=13),axis.text.x=element_text(size=13,angle=angle,vjust=0.5),axis.title=element_text(size=15),legend.text=element_text(size=13),legend.title=element_text(size=15))+
     xlim(score.range[1],score.range[2])+xlab("per-cell CNV variance")
  
  ggsave(paste0(filename,"_density.png"),plot=g1,width=plot_width,height=plot_height,unit="in",device="png")
  ggsave(paste0(filename,"_density.pdf"),plot=g1,width=plot_width,height=plot_height,unit="in",device="pdf")
  
  print("Density plot done!")
  
  g2<-ggplot(dataGG)+
     geom_violin(aes(x=identity,y=score,fill=identity),adjust=violin_adjust)+
     theme(axis.text=element_text(size=13),axis.text.x=element_text(size=13,angle=angle,vjust=0.5),axis.title=element_text(size=15),legend.text=element_text(size=13),legend.title=element_text(size=15))+
     ylab("per-cell CNV variance")
  
  ggsave(paste0(filename,"_violin.png"),plot=g2,width=plot_width,height=plot_height,unit="in",device="png")
  ggsave(paste0(filename,"_violin.pdf"),plot=g2,width=plot_width,height=plot_height,unit="in",device="pdf")
  
  print("Violin plot done!")
  
  return(list(score=score,density_plot=g1,vlnplot=g2))
}



wilcox_test_single<-function(exps,group,ref=1,paired=FALSE,data_type="raw"){
  gr<-levels(group)
  if (ref==1){
    a<-exps[which(group==gr[1])]
    b<-exps[which(group==gr[2])]
  }else{
    a<-exps[which(group==gr[2])]
    b<-exps[which(group==gr[1])]
  }
  if (data_type=="log"){
    res<-mean(b)-mean(a)
  }else{
    res<-mean(b)/mean(a)
  }
  wilcox<-wilcox.test(b,a,paired = paired)
  return(c(res,wilcox$p.value,wilcox$statistic))
}



t_test_single<-function(exps,group,ref=1,paired=FALSE,data_type="raw"){
  gr<-levels(group)
  if (ref==1){
    a<-exps[which(group==gr[1])]
    b<-exps[which(group==gr[2])]
  }else{
    a<-exps[which(group==gr[2])]
    b<-exps[which(group==gr[1])]
  }
  if (data_type=="log"){
    res<-mean(b)-mean(a)
  }else{
    res<-mean(b)/mean(a)
  }
  tt<-t.test(b,a,paired = paired)
  return(c(res,tt$p.value,tt$statistic))
}



group_wise_comparison<-function(compared_value,compared_class_label,test_method="wilcox",test_adjust_method="fdr"){
  ## a convenient function for group pair-wise comparison (especially when there are multiple groups needed to be compared)
  ## currently supported wilcox rank sum test and t-test
  res<-NULL
  classes<-unique(compared_class_label)
  for (i in 1:length(classes)){
    i_index<-which(compared_class_label==classes[i])
    for (j in (i+1):length(classes)){
      if (i==j){
        next
      }
      j_index<-which(compared_class_label==classes[j])
      if (length(i_index)>0 & length(j_index)>0){
        if (test_method=="wilcox"){
          ans<-wilcox.test(compared_value[i_index],compared_value[j_index])
        }
        if (test_method=="t-test"){
          ans<-t.test(compared_value[i_index],compared_value[j_index])
        }
        res_sub<-data.frame(class1=classes[i],class2=classes[j],statistic=ans$statistic,p_val=ans$p.value)
        res<-rbind(res,res_sub)
      }
    }
  }
  res$padj<-p.adjust(res$p_val,test_adjust_method)
  
  res$sig_mark<-""

  sig_mark<-c("*","**","***")
  sig_base<-c(5e-2,1e-2,1e-3)
  names(sig_mark)<-sig_base
  for (i in sig_base){
    res$sig_mark[which(res$padj<i)]<-sig_mark[as.character(i)]
  }
  
  return(res)
}



two_groups_corr<-function(exps1,exps2=NULL,cell_per_experiment=500,experiment=100,random.seed=1709,internal=TRUE,external=TRUE,corr_method="spearman"){
  set.seed(random.seed)
  res_df<-NULL
  if (ncol(exps1)<cell_per_experiment){
    replacement1<-TRUE
  }else{
    replacement1<-FALSE
  }
  
  if (!is.null(exps2)){
    if (ncol(exps2)<cell_per_experiment){
      replacement2<-TRUE
    }else{
      replacement2<-FALSE
    }
  }else{
    external<-FALSE
    internal2_res<-NULL
  }
  
  internal1_res<-numeric(length(experiment))
  internal2_res<-numeric(length(experiment))
  external_res<-numeric(length(experiment))
  
  for (i in 1:experiment){
    z1<-sample(c(1:ncol(exps1)),cell_per_experiment,replace=replacement1)
    x1<-exps1[,z1]
    
    if (!is.null(exps2)){
      z2<-sample(c(1:ncol(exps2)),cell_per_experiment,replace=replacement2)
      x2<-exps2[,z2]
    }
    
    if (internal){
      res<-cor(as.matrix(x1),method=corr_method)
      res<-res[lower.tri(res)]
      internal1_res[i]<-mean(res)
      
      if (!is.null(exps2)){
        res<-cor(as.matrix(x2),method=corr_method)
        res<-res[lower.tri(res)]
        internal2_res[i]<-mean(res)
      }
    }
    
    if (external){
      res<-numeric(length(cell_per_experiment))
      for (j in 1:cell_per_experiment){
        res[j]<-cor(x1[,j],x2[,j],,method=corr_method)
      }
      external_res[i]<-mean(res)
    }
    rm(x1,x2,res)
    gc()
  }
  
  if (!external){
    external_res<-NULL
  }
  if (!internal){
    internal1_res<-NULL
    internal2_res<-NULL
  }
  
  return(list(internal1=internal1_res,internal2=internal2_res,external=external_res))
}



calculate_score_distance<-function(seurat,x.cells,y.cells,cluster_col.x,cluster_col.y,gene_list,
                                   calculate_background=TRUE,aggr_method="median",random.seed=1008,already_calculate_score=FALSE,scale01=TRUE){
  ## gene_list: data.frame include cols "gene", "cluster"
  ## x as query and y as background (reference)
  require(Seurat)
  set.seed(random.seed)
  
  clusters<-unique(gene_list$cluster)
  
  if (!already_calculate_score){
  print("Calculating score...")
    for (cl in clusters){
      print(cl)
      seurat<-AddModuleScore(seurat,features=list(gene_list$gene[which(gene_list$cluster==cl)]),name=cl)
      gc()
    }
  }else{
    print("Score already calculated.")
  }
  
  score_cols<-paste0(clusters,1)
  
  if (scale01){
    meta<-seurat@meta.data
    score_mat<-as.matrix(meta[,score_cols])
    score_mat<-apply(score_mat,2,norm01)
    meta[,score_cols]<-score_mat
    seurat@meta.data<-meta
  }
  
  
  
  x<-subset(seurat,cells=x.cells)
  y<-subset(seurat,cells=y.cells)
  
  meta_x<-x@meta.data
  meta_y<-y@meta.data

  
  if (is.na(cluster_col.x)){  ## per-cell
    xdims<-length(x.cells)
    xcols<-x.cells
    cx_list<-x.cells
  }else{  ## per-group (cluster)
    xcols<-unique(meta_x[,cluster_col.x])
    xdims<-length(xcols)
    cx_list<-xcols
  }
  if (is.na(cluster_col.y)){
    ydims<-length(y.cells)
    ycols<-y.cells
    cy_list<-y.cells
  }else{
    ycols<-unique(meta_y[,cluster_col.y])
    ydims<-length(ycols)
    cy_list<-ycols
  }
  
  if (calculate_background){
    ycols<-c(ycols,"Background")
    ydims<-ydims+1
    cy_list<-ycols
  }
  
  mat<-matrix(nrow=ydims,ncol=xdims)
  colnames(mat)<-xcols
  rownames(mat)<-ycols
  

  print("Aggregate score and calculate distance...")
  for (cy in cy_list){
    # print(paste("Y:",cy))
    if (!(cy %in% clusters) & (cy!="Background")){
      next
    }
    if (is.na(cluster_col.y)){
      sy<-cy
      score_y<-c(as.matrix(meta_y[sy,score_cols]))
      names(score_y)<-score_cols
    }else{
      if (cy=="Background"){
        score_y<-as.matrix(meta_y[,score_cols])
      }else{
        sy<-which(meta_y[,cluster_col.y]==cy)
        score_y<-as.matrix(meta_y[sy,score_cols])
      }
      
      if (aggr_method=="median"){
        score_y<-apply(score_y,2,median)
      }else{
        score_y<-apply(score_y,2,mean)
      }
    }
    
    # print(score_y)
    # print(class(score_y))
    
    if (length(score_y)>1 & cy!="Background"){
      score_y<-score_y[paste0(cy,1)]
    }
    
    for (cx in cx_list){
      # print(paste("X:",cx))
      if (is.na(cluster_col.x)){
        sx<-cx
        score_x<-c(as.matrix(meta_x[sx,score_cols]))
        names(score_x)<-score_cols
      }else{
        sx<-which(meta_x[,cluster_col.x]==cx)
        score_x<-as.matrix(meta_x[sx,score_cols])
        if (aggr_method=="median"){
          score_x<-apply(score_x,2,median)
        }else{
          score_x<-apply(score_x,2,mean)
        }
      }
      
      # print(score_x)
      # print(class(score_x))
      
      if (length(score_x)>1){
        if (cy=="Background"){
          score_x<-score_x[paste0(cx,1)]
          score_y<-score_y[paste0(cx,1)]
        }else{
          score_x<-score_x[paste0(cy,1)]
        }
      }
      mat[cy,cx]<-score_x-score_y

    }
  }
  return(list(mat=mat,seurat=seurat))
}



ssGSEA_wrapped_function<-function(mat,class_label,gmt,permutation_times=100,method="pseudo_bulk",draw_cells=1000,pseudo_bulk_samples=500,
                                  parallel.sz=20,min.sz=10,max.sz=500,ssgsea.norm=TRUE,seed=9953,ssGSEA_mat=NULL){
  require(GSVA)
  require(limma)
  require(Matrix)
  require(parallel)
  
  set.seed(seed)
  
  call_diff<-function(label,mat,lb_order){
    gsva.table<-NULL
    for (cl in lb_order){
      res<-data.frame(cluster=rep("others",ncol(mat)))
      rownames(res)<-colnames(mat)
      res$cluster[which(label==cl)]<-"selected"
      design<-model.matrix(~.,res)

      fit<-lmFit(mat,design)
      fit<-eBayes(fit)
      gsva_res<-topTable(fit,number=Inf,sort.by="none",coef="clusterselected")
      gsva_res$pathway<-rownames(gsva_res)
      gsva_res$cluster<-rep(cl,nrow(gsva_res))
      rownames(gsva_res)<-NULL
      gsva.table<-rbind(gsva.table,gsva_res)
      rm(fit,gsva_res)
      gc()
    }
    return(gsva.table)
  }
  
  generate_pseudo_bulk<-function(mat,class_label,draw_cells,pseudo_bulk_samples,lb_order){
  
    aggr_sum<-function(index,mat){
      return(rowSums(mat[,index]))
    }
  
    bulk_class_label<-NULL
    index.list<-list()
    
    j<-0
    
    for (cl in lb_order){
      bulk_class_label<-c(bulk_class_label,rep(cl,pseudo_bulk_samples))
      zc<-which(class_label==cl)
      if (length(zc)<(draw_cells*1.5)){
        repl<-TRUE
      }else{
        repl<-FALSE
      }
      
      for (i in 1:pseudo_bulk_samples){
        j<-j+1
        index<-zc[sample(length(zc),draw_cells,replace=repl)]
        index.list[[j]]<-index
      }
      
    }
    
    exps<-mclapply(index.list,aggr_sum,mat=mat,mc.cores=parallel.sz,mc.cleanup=TRUE)
    
    pseudo_bulk_mat<-matrix(unlist(exps),nrow=nrow(mat),ncol=(length(lb_order)*pseudo_bulk_samples),byrow=FALSE)
    rownames(pseudo_bulk_mat)<-rownames(mat)
    colnames(pseudo_bulk_mat)<-paste0("p_",c(1:ncol(pseudo_bulk_mat)))
    
    print(dim(pseudo_bulk_mat))
    print(pseudo_bulk_mat[1:5,1:5])
    
    genes_var<-apply(pseudo_bulk_mat,1,var)
    pseudo_bulk_mat<-pseudo_bulk_mat[genes_var>0,]
    
    
    return(list(pseudo_bulk_mat=pseudo_bulk_mat,bulk_class_label=bulk_class_label))
  }
  
  
  
  unique_cl<-unique(class_label)
  kegmt<-GSEABase::getGmt(gmt,geneIdType=GSEABase::SymbolIdentifier())
  
  if (method=="pseudo_bulk"){

    print(Sys.time())
    print("Generating pseudo-bulk expression matrix...")
    
    true_pseudo_bulk<-generate_pseudo_bulk(mat=mat,class_label=class_label,draw_cells=draw_cells,pseudo_bulk_samples=pseudo_bulk_samples,lb_order=unique_cl)
    
    if (is.null(ssGSEA_mat)){
      print(Sys.time())
      print("ssGSEA analysis...")
      ssGSEA_mat<-gsva(true_pseudo_bulk$pseudo_bulk_mat,kegmt,method="ssgsea",parallel.sz=parallel.sz,min.sz=min.sz,max.sz=max.sz,ssgsea.norm=ssgsea.norm)
    }else{
      print("ssGSEA analysis skipped.")
    }
    
    true_res<-call_diff(label=true_pseudo_bulk$bulk_class_label,mat=ssGSEA_mat,lb_order=unique_cl)
    
    true_res$adj.P.Val<-p.adjust(true_res$P.Value)
    
    permutatedLabel.list<-list()
    
    for(j in 1:permutation_times){
      permutatedLabel <- sample(true_pseudo_bulk$bulk_class_label,length(true_pseudo_bulk$bulk_class_label),rep=FALSE) 
      permutatedLabel.list[[j]]<-permutatedLabel
    }
    print(Sys.time())
    print("Calculating FDR by permutation method...")
    
    permuated_res<-mclapply(permutatedLabel.list,call_diff,mat=ssGSEA_mat,lb_order=unique_cl,mc.cores=parallel.sz,mc.cleanup=TRUE)
    
    permuated_B_mat<-matrix(nrow=nrow(true_res),ncol=permutation_times)
    for (j in 1:permutation_times){
      permuated_B_mat[,j]<-permuated_res[[j]]$B
    }
    
    fdr<-numeric(nrow(true_res))
    
    for (j in 1:nrow(true_res)){
      fdr[j]<-length(which(true_res$B[j]<permuated_B_mat[j,]))/ncol(permuated_B_mat)
    }
    
    true_res$FDR<-fdr
    
    return(list(final_res=true_res,permuated_B_mat=permuated_B_mat,ssGSEA_mat=ssGSEA_mat,
                true_pseudo_bulk_exps=true_pseudo_bulk$pseudo_bulk_mat,true_pseudo_bulk_class_label=true_pseudo_bulk$bulk_class_label))
  }
}



microarray_gene_annotation<-function(exp_matrix,cdf.db=NULL,annot.table=NULL,ncores=1){
  ## transfer probeset-level expression to gene-level expression
  ## genes with multiple probeset-level expression values are averaged
  ## return a matrix with genes in rows and samples in cols
  if (!is.null(cdf.db)){
    trans<-select(cdf.db,keys=row.names(exp_matrix),column="SYMBOL",multiVals="asNA")
    trans<-trans[which(!is.na(trans$SYMBOL)),]
    genes<-unique(trans$SYMBOL)
  }else{
    if (!is.null(annot.table)){
      ## annot.table should contain cols "PROBEID" & "SYMBOL"
      trans<-annot.table
      trans<-trans[which(!is.na(trans$SYMBOL) & !(trans$SYMBOL=="")),]
      genes<-unique(trans$SYMBOL)
    }else{
      stop("probeset-gene annotation info not provided, can NOT excute annotation process, exit now...")
    }
  }
  
  avg_exp<-function(gene){
    sel<-trans$PROBEID[which(trans$SYMBOL==gene)]
    if (length(sel)==1){
      exps<-exp_matrix[sel,]
    }else{
      exp_sub<-exp_matrix[sel,]
      exp_sub<-colSums(exp_sub)/nrow(exp_sub)
      exps<-exp_sub
    }
    return(exps)
  }
  if (ncores>1){
    require(parallel)
    exp_matrix_gene<-mclapply(genes,avg_exp, mc.cores = ncores, mc.cleanup = TRUE)
  }else{
    exp_matrix_gene<-lapply(genes,avg_exp)
  }
  
  exp_matrix_gene<-unlist(exp_matrix_gene)
  exp_matrix_gene<-matrix(exp_matrix_gene,ncol=ncol(exp_matrix),byrow=TRUE)
  rownames(exp_matrix_gene)<-genes
  colnames(exp_matrix_gene)<-colnames(exp_matrix)
  return(exp_matrix_gene)
}



ExtractTranscriptLength<-function(gtfFileName,extract_type='exon',group_by='gene_name',generateFile=TRUE){
  require(rtracklayer)

  d<-import(gtfFileName)
  exons<-d[which(d$type==extract_type),]
  tmp <- split(exons,as.character(mcols(exons)[,group_by]))
  Gene_length <- sum(width(reduce(tmp)))
  output<-data.frame(geneSymbol=names(Gene_length),geneLength=Gene_length)
  if (generateFile){
    write.table(output,paste0(gtfFileName,".Length.txt"),sep = "\t",row.names = F,quote = F)
  }
  return(output)
}



volcano_plot_general<-function(df,sig_logfc=1, sig_p_adj=0.05,show_genes=NULL,p_val_cutoff=NA,point_size=1,line_size=1,font_size=2,label_size=0.25){
  ## df: a data.frame with columns: gene, logfc, p, p_adj
  require(ggplot2)
  print(colnames(df))
  
  min_p<-min(df$p_adj[df$p_adj>0])
  df$p_adj[df$p_adj==0]<-min_p*0.1
  df$log_p_adj<-(-log10(df$p_adj))
  
  if (!is.na(p_val_cutoff)){
    print(paste("Using p-value cutoff, p-value smaller than",p_val_cutoff,"will be set to this value"))
    df$log_p_adj[df$log_p_adj>(-log10(p_val_cutoff))]<-(-log10(p_val_cutoff))
  }
  
  df$color<-"grey"
  df$color[which(df$logfc>sig_logfc & df$p_adj<sig_p_adj)]<-"red"
  df$color[which(df$logfc<(-sig_logfc) & df$p_adj<sig_p_adj)]<-"navyblue"
  
  df$label<-""
  df$label[df$gene %in% show_genes]<-df$gene[df$gene %in% show_genes]
  
  g<-ggplot(df)+
     geom_point(aes(x=logfc,y=log_p_adj,color=color),size=point_size)+
	 scale_color_manual(values=c("red"="red","navyblue"="navyblue","grey"="grey"))+
	 geom_vline(xintercept=sig_logfc,linetype='dashed',color="#595959",size=line_size)+
	 geom_vline(xintercept=-sig_logfc,linetype='dashed',color="#595959",size=line_size)+
	 geom_hline(yintercept=(-log10(sig_p_adj)),linetype='dashed',color="#595959",size=line_size)+
	 xlab("log Fold Change of gene expression")+ylab('-log10 adjusted P-value')
	 
  if (!is.null(show_genes)){
    require(ggrepel)
	g<-g+geom_text_repel(aes(x=logfc,y=log_p_adj,label=label),max.overlaps=Inf,min.segment.length=0,size=font_size,direction='both', label.size = label_size,max.time=5,max.iter=1e7)
  }
  
  g<-g+theme_classic()+theme(legend.position='none')
  
  return(g)
}



plot_pseudotime_heatmap_CH<-function(exps,meta,binned_col='Pseudotime',bins=30,data_min=-Inf,data_max=Inf,color=NULL,
                                     cluster_rows=FALSE,scale_method='none',centering='none',show_row_names=TRUE,row_dend_width =  unit(10, "mm"),
									 dist_method='euclidean',hclust_method='ward.D2',cluster_num=1,random_seed=183,
									 heatmap_legend_param=NULL,
									 smoothed=FALSE,smooth_span=0.7){
  require(ComplexHeatmap)
  require(dendextend)
  
  set.seed(random_seed)
  
  ## this is a mimic function of 'plot_pseudotime_heatmap' in 'monocle' utilizing ComplexHeatmap
  ## exps and meta are assumed to be already aligned (also, ncol of exps and nrow of meta should be the same)
  
  ## centering -> scaling -> trimming -> smoothing
  
  meta$cut_bins<-cut(meta[,binned_col],breaks=bins)
  
  new_mat<-matrix(nrow=nrow(exps),ncol=bins)
  colnames(new_mat)<-levels(meta$cut_bins)
  rownames(new_mat)<-rownames(exps)
  for (i in levels(meta$cut_bins)){
    exps_sub<-exps[,meta$cut_bins==i]
	if (is.null(dim(exps_sub))){
	  RMeans<-exps_sub
	}else{
	  RMeans<-rowMeans(exps_sub,na.rm=TRUE)
	}
	new_mat[,i]<-RMeans
  }
  
  if (scale_method=='row'){
    new_mat<-t(scale(t(new_mat)))
  }else{
    if (scale_method=='col'){
	  new_mat<-scale(new_mat)
	}else{
	  if (scale_method=='both'){
	    ## first row, then col
	    new_mat<-t(scale(t(new_mat)))
		new_mat<-scale(new_mat)
	  }
	}
  }
  
  if (scale_method!='none'){
    centering<-'none'
  }
  if (centering=='row'){
    new_mat<-t(scale(t(new_mat),scale=FALSE))
  }else{
    if (centering=='col'){
	  new_mat<-scale(new_mat,scale=FALSE)
	}else{
	  if (centering=='both'){
	    ## first row, then col
	    new_mat<-t(scale(t(new_mat),scale=FALSE))
		new_mat<-scale(new_mat,scale=FALSE)
	  }
	}
  }
  
  if (is.finite(data_min)){
    new_mat[new_mat<data_min]<-data_min
  }
  if (is.finite(data_max)){
    new_mat[new_mat>data_max]<-data_max
  }
  
  new_mat_orig<-new_mat
  
  if (smoothed){
    packed_func<-function(y,span){
	  d<-data.frame(x=c(1:length(y)),y=y)
	  loess_res <- loess(y ~ x, data=d, span=span)
      pred <- predict(loess_res)
	  return(pred)
	}
	
    smoothed_mat<-apply(new_mat,1,packed_func,span=smooth_span)
	smoothed_mat<-t(smoothed_mat)
	new_mat<-smoothed_mat
	colnames(new_mat)<-colnames(new_mat_orig)
	rownames(new_mat)<-rownames(new_mat_orig)
  }else{
    smoothed_mat<-NULL
  }
    
  
  if (cluster_rows){
    row_dend = as.dendrogram(hclust(dist(new_mat,method=dist_method),method=hclust_method))
    row_dend = color_branches(row_dend, k = cluster_num) # `color_branches()` returns a dendrogram object
  }else{
    row_dend<-FALSE
  }
  
  if (is.null(heatmap_legend_param)){
    heatmap_legend_param<-list(title='Expression',grid_height=unit(8,'mm'),title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6))
  }
  
  if (is.null(color)){
    color = c("blue", "#EEEEEE", "red")
  }
  
  Hob<-Heatmap(new_mat,cluster_rows=row_dend,col=color,
               row_names_side = 'right', row_names_gp = gpar(fontsize = 6),show_row_names=show_row_names,
               cluster_columns=FALSE,column_order=levels(meta$cut_bins),show_column_names =FALSE,
			   heatmap_legend_param = heatmap_legend_param)
  return(list(g=Hob,mat=new_mat_orig,smoothed_mat=smoothed_mat))
}
