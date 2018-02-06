## Exercise 01 - basics of Snakemake rules

**Goal:** turn a shell script (`aln.sh` - alignment with hisat2) into a single rule
Snakemake file (`Snakefile`). One possible solution is given in `Snakefile.finished`.


The script `aln.sh` takes a fastq file as an argument and aligns it to the
S. cerevisiae genome:

```ShellSession
[user@cn1234]$ cat aln.sh

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

[user@cn1234]$ ./aln.sh 00fastq/ERR458495.fastq.gz

1066501 reads; of these:
  1066501 (100.00%) were unpaired; of these:
    52263 (4.90%) aligned 0 times
    786836 (73.78%) aligned exactly 1 time
    227402 (21.32%) aligned >1 times
95.10% overall alignment rate
```


