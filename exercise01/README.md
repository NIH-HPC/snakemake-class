## Exercise 01 - basics of Snakemake rules

**Goal:** turn a shell script (`aln.sh` - alignment with hisat2) into a single
rule Snakemake file (`Snakefile`). One possible solution is given in
**Snakefile.finished**.


:information_source: Snakefiles are basically Python with some extra constructs for defining
rules based workflows. That means that **Snakefiles are whitesspace sensitive** - 
indentation matters and tabs and spaces can't be mixed. Code in these examples
uses spaces - no tabs. Please adjust your editors accordingly.

The script `aln.sh` takes a fastq file as an argument and aligns it to the
S. cerevisiae genome. We use the tools installed in the container that
was downloaded during the setup. If you haven't done so in the previous
exercise please set up all the required bind mounts with

```console
user@cn1234> source /usr/local/current/singularity/app_conf/sing_binds
```

Then run the script for one sample

```console
user@cn1234> cat aln.sh
#! /bin/bash

fq=$1
bam=02aln/$(basename $fq .fastq.gz).bam 

idx=00ref/hisat_index/R64-1-1
mkdir -p 02aln
hisat2 -k 4 -x $idx -U $fq --threads 4 \
  | samtools sort -T $bam -O BAM \
  > $bam
samtools index $bam

user@cn1234> export CONTAINER="$(cd .. && pwd)/00container/rnaseq.sif"
user@cn1234> singularity exec "$CONTAINER" ./aln.sh 00fastq/ERR458495.fastq.gz

1066501 reads; of these:
  1066501 (100.00%) were unpaired; of these:
    37947 (3.56%) aligned 0 times
    912472 (85.56%) aligned exactly 1 time
    116082 (10.88%) aligned >1 times
96.44% overall alignment rate

user@cn1234> rm -rf 02aln
```

Let's turn this into a simple Snakefile by formalizing the input and
output file types. Note that by convention snakemake looks for a file named
`Snakefile` if it is not given a different filename with the `-s` option.

A simple snakefile with a general rule for how to create bam files in the
`02aln` directory from compressed fastq files in the `00fastq` directory might
look similar to the following:

```python
# this file is 'Snakefile'
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    shell:
        """
        hisat2 -k 4 -x {input.idx} -U {input.fq} --threads 4 \
          | samtools sort -T {output.bam} -O BAM \
          > {output.bam}
        samtools index {output.bam}
        """
```

When snakemake is asked to create a specific file it searches for a rule with
an output pattern that matches the file name and then executes the rule, which
in this example is shell code. The input and output file names for the specific
file are available as string substitutions within the shell code section. Note that
specifying `--cores` is required.

```console
user@cn1234> snakemake --cores=4 -np 02aln/ERR458502.bam
Building DAG of jobs...
Job stats:
job       count
------  -------
hisat2        1
total         1


[Thu Sep 26 17:13:36 2024]
rule hisat2:
    input: 00fastq/ERR458502.fastq.gz, 00ref/hisat_index/R64-1-1
    output: 02aln/ERR458502.bam, 02aln/ERR458502.bam.bai
    jobid: 0
    reason: Missing output files: 02aln/ERR458502.bam
    wildcards: sample=ERR458502
    resources: tmpdir=/tmp


        hisat2 -k 4 -x 00ref/hisat_index/R64-1-1 -U 00fastq/ERR458502.fastq.gz --threads 4           | samtools sort -T tmp/ERR458502 -O BAM           > 02aln/ERR458502.bam
        samtools index 02aln/ERR458502.bam

Job stats:
job       count
------  -------
hisat2        1
total         1

Reasons:
    (check individual jobs above for details)
    missing output files:
        hisat2

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

This can't actually be run yet because hisat and samtools are not
available on our path. However we can teach snakemake to run a particular
rule inside a singularity container by making one small change:

```python
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    singularity:
        "library://wresch/classes/rnaseq:0.8"
    shell:
        """
        hisat2 -k 4 -x {input.idx} -U {input.fq} --threads 4 \
          | samtools sort -T tmp/{wildcards.sample} -O BAM \
          > {output.bam}
        samtools index {output.bam}
        """
```

By default, snakemake ignores the `singularity` directive. The
`--use-singularity` option is required to enable use of singularity and the
`singularity` executable has to be available on the path. Additional
singularity options can be passed with `--singularity-args` and the
location of pulled containers can be specified with `--singularity-prefix` which
defaults to `.snakemake/singularity`. We already have the container pulled
down in `../00container/` so we'll use that instead of creating another copy.
We are using the environment variable `$SINGULARITY_BINDPATH` so
for now we don't need `--singularity-args`


```console
user@cn1234> snakemake --cores 8 --use-singularity 02aln/ERR458502.bam \
    --singularity-prefix=../00container

Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job stats:
job       count
------  -------
hisat2        1
total         1

Select jobs to execute...

[Thu Sep 26 17:19:17 2024]
rule hisat2:
    input: 00fastq/ERR458502.fastq.gz, 00ref/hisat_index/R64-1-1
    output: 02aln/ERR458502.bam, 02aln/ERR458502.bam.bai
    jobid: 0
    reason: Missing output files: 02aln/ERR458502.bam
    wildcards: sample=ERR458502
    resources: tmpdir=/tmp

Activating singularity image /spin1/users/wresch/code/class_materials/snakemake-class/00container/2354d2ff28bcf0b42c57fae398b4c9b5.simg
1853031 reads; of these:
  1853031 (100.00%) were unpaired; of these:
    49450 (2.67%) aligned 0 times
    1655519 (89.34%) aligned exactly 1 time
    148062 (7.99%) aligned >1 times
97.33% overall alignment rate
[Thu Sep 26 17:19:29 2024]
Finished job 0.
1 of 1 steps (100%) done
Complete log: .snakemake/log/2024-09-26T171916.905065.snakemake.log


```

If we ask snakemake to produce the same output file again, it recognizes
that nothing has to be done because the output files are newer than
the input files:

```console
user@cn1234> snakemake --cores 8 --use-singularity 02aln/ERR458502.bam \
    --singularity-prefix=../00container
Building DAG of jobs...
Nothing to be done.
```

If the modification time on one of the input files is changed to a more
recent time than the output files with `touch`, the output files are
regenerated. `snakemake --summary` can be used to get a summary of 
all outputs and their status.

```console
user@cn1234> touch 00fastq/ERR458502.fastq.gz
user@cn1234> snakemake --summary 02aln/ERR458502.bam
Building DAG of jobs...
output_file     date    rule    version log-file(s)     status  plan
02aln/ERR458502.bam     Thu Sep 26 17:19:29 2024        hisat2  -               updated input files     update pending
02aln/ERR458502.bam.bai Thu Sep 26 17:19:29 2024        hisat2  -               updated input files     update pending


user@cn1234> snakemake --cores 8 --use-singularity 02aln/ERR458502.bam \
    --singularity-prefix=../00container
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

Since snakemake 7.8 reruns are also triggered if parameters, code, input file
set, or software stack changed. For example if we make a small change in the 
hisat2 rule by adding a new empty line we can see that

```console
user@cn1234> snakemake --summary 02aln/ERR458502.bam
Building DAG of jobs...
output_file     date    rule    version log-file(s)     status  plan
02aln/ERR458502.bam     Thu Sep 26 17:21:58 2024        hisat2  -               rule implementation changed     update pending
02aln/ERR458502.bam.bai Thu Sep 26 17:21:58 2024        hisat2  -               rule implementation changed     update pending
```


Now let's add a rule for cleaning up all generated files at the end of the
snakefile.  A rule to clean up generated files should in general not be the
first rule in a Snakefile since the first rule is the default rule snakemake
runs when no rule or concrete file is provided as an argument.

```python
rule clean:
    shell:
        """
        rm -rf 02aln
        """
```

Snakemake can be asked to run a rule as long as it does not have any wildcard
inputs:

```console
user@cn1234> snakemake -pn clean
Building DAG of jobs...

rule clean:
    jobid: 0


        rm -rf 02aln
```

This does not work with a rule that has wildcard inputs b/c that is an abstract
rule and snakemake would have no way to know what the actual input and output
files should be:

```console
user@cn1234> snakemake -pn hisat2
WorkflowError:                                                                                    
Target rules may not contain wildcards. Please specify concrete files or a rule without wildcards.
```

So, we have to call snakemake with multiple targets to generate multiple
alignments:

```console
user@cn1234> snakemake --cores 8 --use-singularity \
    --singularity-prefix=../00container \
    02aln/ERR458495.bam 02aln/ERR458502.bam
```

Or better - create a default rule (i.e. first rule in the file) that lists
each of the desired output files as an input but does not do any actual work.
Such a rule is conventionally called `all`:

```console
rule all:
    input: "02aln/ERR458495.bam",
           "02aln/ERR458502.bam"
```

Now, if we do a `snakemake all`, or just `snakemake`, snakemake will run the
all rule which requires two inputs but has no action. It then searches for
rules to generate those input files and, in this case, finds that hisat2 can generate
both of the input files and executes the rule twice with the different input files:

```console
user@cn1234> snakemake --use-singularity --cores=8 --singularity-prefix=../00container all
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       all
        2       hisat2
        3
...
```


