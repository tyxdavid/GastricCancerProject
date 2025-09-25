source('./Utils.R')

d<-readRDS("exps.rds")

anno_table<-read.table('./official_symbol.txt',sep="\t",header=TRUE)
anno_table$SYMBOL<-anno_table$SYMBOL

exp_matrix_gene<-microarray_gene_annotation(exp_matrix=d,annot.table=anno_table,ncores=20)

saveRDS(exp_matrix_gene,'exp_matrix_gene.rds')
