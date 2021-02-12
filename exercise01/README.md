## Exercise 01 - basics of Snakemake rules

**Goal:** turn a shell script (`aln.sh` - alignment with hisat2) into a single
rule Snakemake file (`Snakefile`). One possible solution is given in
**Snakefile.finished**.


:information_source: Snakefiles are basically Python with some extra constructs for defining
rules based workflows. That means that **Snakefiles are whitesspace sensitive** - 
indentation matters and tabs and spaces can't be mixed. Code in these examples
uses spaces - no tabs. Please adjust your editors accordingly.

The script `aln.sh` takes a fastq file as an argument and aligns it to the
S. cerevisiae genome. We use the tools installed in the container via the
rnaseq wrapper introduced in the previous example

```console
user@cn1234> cat aln.sh

fq=$1
bam=02aln/$(basename $fq .fastq.gz).bam 

idx=00ref/hisat_index/R64-1-1
mkdir -p 02aln
hisat2 -k 4 -x $idx -U $fq --threads 4 \
  | samtools sort -T $bam -O BAM \
  > $bam
samtools index $bam

user@cn1234> ./rnaseq aln.sh 00fastq/ERR458495.fastq.gz

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
file are available as string substitutions within the shell code section.

```console
user@cn1234> snakemake -np 02aln/ERR458502.bam

rule hisat2:
    input: 00fastq/ERR458502.fastq.gz, 00ref/hisat_index/R64-1-1
    output: 02aln/ERR458502.bam, 02aln/ERR458502.bam.bai
    jobid: 0
    wildcards: sample=ERR458502


        hisat2 -k 4 -x 00ref/hisat_index/R64-1-1 -U 00fastq/ERR458502.fastq.gz --threads 4 \
        | samtools sort - T 02aln/ERR458502.bam -O BAM \
        > 02aln/ERR458502.bam
        samtools index 02aln/ERR458502.bam
```

This can't actually be run yet because hisat and samtools are not
available on the path. However we can teach snakemake to run a particular
rule inside a singularity container by making one small change:

```python
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    singularity:
        "library://wresch/classes/rnaseq:0.5"
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
singularity options can be passed with `--singularity-args`. We will use this
to pass in bind mount options. Finally, we want to store containers in the
`00container` directory in the root dir of the repo using the
`--singularity-prefix` option:

```console
user@cn1234> snakemake --cores 8 --use-singularity 02aln/ERR458502.bam

Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 8
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       hisat2
        1

[Fri Feb 12 10:50:27 2021]
rule hisat2:
    input: 00fastq/ERR458502.fastq.gz, 00ref/hisat_index/R64-1-1
    output: 02aln/ERR458502.bam, 02aln/ERR458502.bam.bai
    jobid: 0
    wildcards: sample=ERR458502

Activating singularity image /spin1/users/user/class_materials/snakemake-class/exercise01/.snakemake/singularity/b891a3a5ef0d77720e876a3540e08ee7.simg
1853031 reads; of these:
  1853031 (100.00%) were unpaired; of these:
    49450 (2.67%) aligned 0 times
    1655519 (89.34%) aligned exactly 1 time
    148062 (7.99%) aligned >1 times
97.33% overall alignment rate
[Fri Feb 12 10:50:49 2021]
Finished job 0.
1 of 1 steps (100%) done
Complete log: /spin1/users/user/class_materials/snakemake-class/exercise01/.snakemake/log/2021-02-12T105027.048281.snakemake.log

```

The singularity container will be cached in the .snakemake subdirectory
of the working directory.

If we ask snakemake to produce the same output file again, it recognizes
that nothing has to be done because the output files are newer than
the input files:

```console
user@cn1234> snakemake --cores 8 --use-singularity 02aln/ERR458502.bam
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
output_file     date    rule    version log-file(s)     status  plan
02aln/ERR458502.bam     Fri Feb 12 10:50:49 2021        hisat2  -               updated input files   update pending
02aln/ERR458502.bam.bai Fri Feb 12 10:50:49 2021        hisat2  -               updated input files   update pending


user@cn1234> snakemake --cores 8 --use-singularity 02aln/ERR458502.bam
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
user@cn1234> snakemake --use-singularity --cores=8 all
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


