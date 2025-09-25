#!/bin/sh

for k in $(cat $1)
do
{
## retain unique mapped reads and sort
samtools view -@ 12 -q 30 -b ./aligned/${k}.unsorted.sam | samtools sort -@ 54 -o ./aligned/${k}.sorted.bam
samtools index -@ 8 ./aligned/${k}.sorted.bam

## remove reads from chrM and scaffolds
samtools view -@ 12 -q 30 -O SAM -h ./aligned/${k}.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort -@ 52 -o ./aligned/${k}.sorted.rmChrM.bam
samtools index -@ 8 ./aligned/${k}.sorted.rmChrM.bam

## mark and remove duplicated reads
java -XX:ParallelGCThreads=32 -Djava.io.tmpdir=/tmp -jar /data1/software/picard.jar MarkDuplicates \
  QUIET=false INPUT=./aligned/${k}.sorted.rmChrM.bam OUTPUT=./aligned/${k}.sorted.rmDup.bam METRICS_FILE=./aligned/${k}.rmDup.metrics \
  REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp

}
done
