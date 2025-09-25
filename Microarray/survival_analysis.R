source("./Utils.R")

library(ggplot2)
library(survival)
library(survminer)
library(survcomp)


exps<-readRDS('./exp_matrix_gene.rds')
meta<-readRDS('./sample_meta.rds')

exps<-exps[,colnames(exps) %in% meta$GEO_ID]
z<-factor(colnames(exps),levels=meta$GEO_ID,ordered = T)
exps<-exps[,order(z)]


meta$group<-'Low'
meta$group[meta$gene_signature_score>median(meta$gene_signature_score)]<-'High'

meta$CEBPB_exps<-exps['CEBPB',]
saveRDS(meta,'meta_v3.rds')


fit <- survfit(Surv(time=OS_month, event=OS) ~ group, data = meta)
ggsurv <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = meta,             # data used to fit survival curves.
  risk.table = F,       # show risk table.
  pval = TRUE, # show p-value of log-rank test.
  pval.size=2,
  conf.int = FALSE,         # show confidence intervals for 
  censor.shape=124,
  censor.size=1,
  palette = c("#BB0021", "#3B4992"),
  xlab = "Time (Months)",   # customize X axis label.
  break.time.by = 24,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(base_size = 6)+theme(axis.title=element_text(size=6)), # customize plot and risk table with a theme.
  tables.theme=theme_classic(base_size = 6)+theme(axis.title=element_text(size=6)),
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  risk.table.fontsize=2,
  linewidth=0.005,
  conf.int.style = "step",  # customize style of confidence intervals
  legend='none',
  newpage=F
)

pdf('GSig_OS.pdf', width = 2.1654, height = 2.0866)
print(ggsurv)
dev.off()


