#!/bin/sh


macs2 callpeak -t ../aligned/NCI-N87-1.sorted.shifted.bam \
-f BAMPE --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all \
--outdir ./NCI-N87-1 -n NCI-N87-1

macs2 callpeak -t ../aligned/NCI-N87-2.sorted.shifted.bam \
-f BAMPE --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all \
--outdir ./NCI-N87-2 -n NCI-N87-2

