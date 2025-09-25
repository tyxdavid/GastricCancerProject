#!/bin/sh

trim_galore --fastqc --fastqc_args " --outdir ./fastqc/" \
-o ./trimmed --cores 4 --paired ./fastq/NCI-N87-1.Rawdata.R1.fq.gz ./fastq/NCI-N87-1.Rawdata.R2.fq.gz

trim_galore --fastqc --fastqc_args " --outdir ./fastqc/" \
-o ./trimmed --cores 4 --paired ./fastq/NCI-N87-2.Rawdata.R1.fq.gz ./fastq/NCI-N87-2.Rawdata.R2.fq.gz
