#! /bin/bash

fq=$1
bam=02aln/$(basename $fq .fastq.gz).bam 

module load hisat/2.0.5 samtools/1.4
idx=00ref/hisat_index/R64-1-1
mkdir -p 02aln
hisat2 -k 4 -x $idx -U $fq --threads 4 \
  | samtools sort -T tmp/ERR458495 -O BAM \
  > $bam
samtools index $bam
