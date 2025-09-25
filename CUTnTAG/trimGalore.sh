#!/bin/sh

trim_galore --fastqc --fastqc_args " --outdir ./fastqc/" \
-o ./trimmed --cores 4 --paired ./fastq/NCI-N87-CEBPB-1_R1.fastq.gz ./fastq/NCI-N87-CEBPB-1_R2.fastq.gz

trim_galore --fastqc --fastqc_args " --outdir ./fastqc/" \
-o ./trimmed --cores 4 --paired ./fastq/NCI-N87-CEBPB-2_R1.fastq.gz ./fastq/NCI-N87-CEBPB-2_R2.fastq.gz
