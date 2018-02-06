## Exercise 01 - basics of Snakemake rules

**Goal:** turn a shell script (`aln.sh` - alignment with hisat2) into a single rule
Snakemake file (`Snakefile`). One possible solution is given in `Snakefile.finished`.


The script `aln.sh` takes a fastq file as an argument and aligns it to the
S. cerevisiae genome:

```console
user@cn1234> cat aln.sh

#! /bin/bash

fq=$1
bam=02aln/$(basename $fq .fastq.gz).bam

module load hisat/2.0.5 samtools/1.4
idx=00ref/hisat_index/R64-1-1
mkdir -p 02aln
hisat2 -k 4 -x $idx -U $fq --threads 4 \
  | samtools sort -T $bam -O BAM \
  > $bam
samtools index $bam

user@cn1234> bash aln.sh 00fastq/ERR458495.fastq.gz

1066501 reads; of these:
  1066501 (100.00%) were unpaired; of these:
    52263 (4.90%) aligned 0 times
    786836 (73.78%) aligned exactly 1 time
    227402 (21.32%) aligned >1 times
95.10% overall alignment rate
```

Let's turn this into a simple Snakefile by formalizing the input and
output. Note that by convention snakemake looks for a file named
`Snakefile` if it is not given a different filename with the `-s` option.

The snakefile contains a single rule which describes how to create a specific
type of output file (or multiple) from a specific type of input file. Here the
inputs are a fastq file and a hisat2 index and the output a bam file and an
index.

```python
# this file is 'Snakefile'
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    shell:
        """
        module load hisat/2.0.5 samtools/1.4
        hisat2 -k 4 -x {input.idx} -U {input.fq} --threads 4 \
          | samtools sort -T {output.bam} -O BAM \
          > {output.bam}
        samtools index {output.bam}
        """
```

When snakemake is asked to create a concrete file it searches for a rule which
has an output that matches this concrete file and then executes the shell code.

```console
user@cn1234> snakemake 02aln/ERR458502.bam
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       hisat2
        1

rule hisat2:
    input: 00fastq/ERR458502.fastq.gz, 00ref/hisat_index/R64-1-1
    output: 02aln/ERR458502.bam, 02aln/ERR458502.bam.bai
    jobid: 0
    wildcards: sample=ERR458502

[+] Loading hisat 2.0.5 on cn3140
[+] Loading samtools 1.4 ...
1853031 reads; of these:
  1853031 (100.00%) were unpaired; of these:
    67169 (3.62%) aligned 0 times
    1428567 (77.09%) aligned exactly 1 time
    357295 (19.28%) aligned >1 times
96.38% overall alignment rate
Finished job 0.
1 of 1 steps (100%) done
```

If we ask snakemake to produce the same output file again, it recognizes
that nothing has to be done because the output files are newer than
the input files:

```console
user@cn1234> snakemake 02aln/ERR458502.bam
Building DAG of jobs...
Nothing to be done.
```

If that changes, the output files are regenerated:
```console
user@cn1234> touch 00fastq/ERR458502.fastq.gz
user@cn1234> snakemake --summary 02aln/ERR458502.bam
output_file     date    rule    version log-file(s)     status  plan                       02aln/ERR458502.bam     Mon Feb  5 22:12:02 2018        hisat2  -               updated inputfiles   update pending
02aln/ERR458502.bam.bai Mon Feb  5 22:12:02 2018        hisat2  -               updated inputfiles   update pending

user@cn1234> snakemake 02aln/ERR458502.bam
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       hisat2
        1
...
```

