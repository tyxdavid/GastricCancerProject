#!/bin/sh


bowtie2 -p 42 \
--very-sensitive --no-discordant -X 2000 \
-x ./genome/bowtie2/genome \
-1 ./trimmed/NCI-N87-1.Rawdata.R1_val_1.fq.gz \
-2 ./trimmed/NCI-N87-1.Rawdata.R2_val_2.fq.gz \
-S ./aligned/NCI-N87-1.unsorted.sam

bowtie2 -p 42 \
--very-sensitive --no-discordant -X 2000 \
-x ./genome/bowtie2/genome \
-1 ./trimmed/NCI-N87-2.Rawdata.R1_val_1.fq.gz \
-2 ./trimmed/NCI-N87-2.Rawdata.R2_val_2.fq.gz \
-S ./aligned/NCI-N87-2.unsorted.sam
