#!/bin/sh


LBpath="./CUT-TAG/CEBPB/NCI-N87/aligned"

sample="NCI-N87-CEBPB-1.sorted.rmDup"
sh SEACR_1.3.sh ${LBpath}/${sample}.fragments.bedgraph ./CUT-TAG/IgG/NCI-N87/aligned/NCI-N87-merged-IgG.sorted.rmDup.fragments.bedgraph non stringent NCI-N87-CEBPB-1

sample="NCI-N87-CEBPB-2.sorted.rmDup"
sh SEACR_1.3.sh ${LBpath}/${sample}.fragments.bedgraph ./CUT-TAG/IgG/NCI-N87/aligned/NCI-N87-merged-IgG.sorted.rmDup.fragments.bedgraph non stringent NCI-N87-CEBPB-2
