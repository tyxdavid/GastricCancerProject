#!/bin/sh

for k in $(cat $1)
do
{
bowtie2 -p 54 \
--very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
-x /data1/liangjq/GC_0622/genome/bowtie2/genome \
-1 ./trimmed/${k}_R1_val_1.fq.gz \
-2 ./trimmed/${k}_R2_val_2.fq.gz \
-S ./aligned/${k}.unsorted.sam

}
done
