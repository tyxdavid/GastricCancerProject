library(ATACseqQC)
setwd('/data1/liangjq/ATAC-seq/aligned')

tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
# library(BSgenome.Hsapiens.UCSC.hg38)
# which <- as(seqinfo(Hsapiens), "GRanges")

bamfile <- "NCI-N87-1.sorted.rmBlack.bam"
gal <- readBamFile(bamfile, tag=tags, asMates=TRUE, bigFile=FALSE)
objs <- shiftGAlignmentsList(gal,outbam="NCI-N87-1.sorted.shifted.bam")

bamfile <-  "NCI-N87-2.sorted.rmBlack.bam"
gal <- readBamFile(bamfile, tag=tags, asMates=TRUE, bigFile=FALSE)
objs <- shiftGAlignmentsList(gal,outbam="NCI-N87-2.sorted.shifted.bam")
