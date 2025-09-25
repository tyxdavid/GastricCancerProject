#!/bin/sh

samtools view -@ 8 -c NCI-N87-CEBPB-1.sorted.rmChrM.bam
samtools view -@ 8 -c NCI-N87-CEBPB-2.sorted.rmChrM.bam

bedtools genomecov -bg -scale 0.017743 -ibam NCI-N87-CEBPB-1.sorted.rmChrM.bam >NCI-N87-CEBPB-1.sorted.rmChrM.bedgraph &
bedtools genomecov -bg -scale 0.034097 -ibam NCI-N87-CEBPB-2.sorted.rmChrM.bam >NCI-N87-CEBPB-2.sorted.rmChrM.bedgraph &
