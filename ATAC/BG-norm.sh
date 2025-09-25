#!/bin/sh

bedtools multicov -p -bams NCI-N87-1.sorted.shifted.bam -bed ./genome/genes/genes_TSS_filtered.bed >NCI-N87-1.sorted.shifted.TSS.count.tsv
bedtools multicov -p -bams NCI-N87-2.sorted.shifted.bam -bed ./genome/genes/genes_TSS_filtered.bed >NCI-N87-2.sorted.shifted.TSS.count.tsv

bedtools genomecov -bg -scale 0.142857 -ibam NCI-N87-1.sorted.shifted.bam >NCI-N87-1.sorted.shifted.bedgraph
bedtools genomecov -bg -scale 0.124515 -ibam NCI-N87-2.sorted.shifted.bam >NCI-N87-2.sorted.shifted.bedgraph

