#!/bin/sh

samtools view -@ 12 -q 30 -O SAM -h ./NCI-N87/aligned/NCI-N87-1.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort -@ 42 -o ./NCI-N87/aligned/NCI-N87-1.sorted.rmChrM.bam
samtools view -@ 12 -q 30 -O SAM -h ./NCI-N87/aligned/NCI-N87-2.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort -@ 42 -o ./NCI-N87/aligned/NCI-N87-2.sorted.rmChrM.bam


java -XX:ParallelGCThreads=32 -Djava.io.tmpdir=/tmp -jar /data1/software/picard.jar MarkDuplicates \
  QUIET=false INPUT=./NCI-N87/aligned/NCI-N87-1.sorted.rmChrM.bam OUTPUT=./NCI-N87/aligned/NCI-N87-1.sorted.rmDup.bam METRICS_FILE=./NCI-N87/aligned/NCI-N87-1.rmDup.metrics \
  REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp

java -XX:ParallelGCThreads=32 -Djava.io.tmpdir=/tmp -jar /data1/software/picard.jar MarkDuplicates \
  QUIET=false INPUT=./NCI-N87/aligned/NCI-N87-2.sorted.rmChrM.bam OUTPUT=./NCI-N87/aligned/NCI-N87-2.sorted.rmDup.bam METRICS_FILE=./NCI-N87/aligned/NCI-N87-2.rmDup.metrics \
  REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp

samtools view -@ 8 -b -L /data1/public_data/hg38-blacklist.v2.bed -U ./NCI-N87/aligned/NCI-N87-1.sorted.rmBlack.bam ./NCI-N87/aligned/NCI-N87-1.sorted.rmDup.bam > overlappingSpecificRegions.bam
samtools index -@ 12 ./NCI-N87/aligned/NCI-N87-1.sorted.rmBlack.bam
samtools view -@ 8 -b -L /data1/public_data/hg38-blacklist.v2.bed -U ./NCI-N87/aligned/NCI-N87-2.sorted.rmBlack.bam ./NCI-N87/aligned/NCI-N87-2.sorted.rmDup.bam > overlappingSpecificRegions-2.bam
samtools index -@ 12 ./NCI-N87/aligned/NCI-N87-2.sorted.rmBlack.bam
