## Exercise 02 - Parallelizing pipelines (1)

**Goal:** Formalize the number of CPUs each rule requires for each rule as well as 
parameters for tools used. Run workflow in parallel. Use the file `Snakefile` as a starting
point. `Snakefile.finished` is a possible solution.

### Threads

In order to run individual tasks in parallel, snakemake needs to know how many
CPUs each rule would like to use and how many CPUs are available in total. In
our example, we would like to run hisat2 with 4 threads. If we allow snakemake
to use 8 CPUs, it could run two alignments in parallel.  If fewer than 4 CPUs
were available, the number of threads for the rule would be scaled down
accordingly. The number of CPUs each rule would like to use is specified
in the `threads` section and is available in the shell block as `{threads}`.
For example, here is the modified hisat2 rule:


```python
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    threads: 4
    singularity:
        "library://wresch/classes/rnaseq:0.6"
    shell:
        """
        mkdir -p 02aln
        hisat2 -k 4 -x {input.idx} -U {input.fq} --threads {threads} \
          | samtools sort -T tmp/{wildcards.sample} -O BAM \
          > {output.bam}
        samtools index {output.bam}
        """
```

Now, the workflow can run in parallel:

```console
user@cn1234> snakemake --cores 8 --use-singularity --singularity-prefix=../00container \
    --singularity-args '-B $PWD:/data --pwd /data'
```

### Resources

In the resources section arbitrary resources (other than threads) required by
the rule can be specified. This can be used to specify, for example, the amount
of memory used by the rule (and in fact the resource name `mem_mb` is special
b/c it is used by the kubernetes runner). But it can also be used to limit the
number of concurrent I/O intensive jobs or include the walltime limit. The total
amount of resources available can be set on the command line. All resources must be
integers.

```python
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    threads: 4
    resources: mem_mb = 6144
    singularity:
        "library://wresch/classes/rnaseq:0.6"
    shell:
        """
        mkdir -p 02aln
        hisat2 -k 4 -x {input.idx} -U {input.fq} --threads {threads} \
          | samtools sort -T tmp/{wildcards.sample} -O BAM \
          > {output.bam}
        samtools index {output.bam}
        """
```

Note that there are different ways in which resources such as memory can
be specified. For example cluster-specific resources could also be specified
in a separate config file which will be discussed later.

```console
user@cn1234> snakemake --cores 8 --resources mem_mb=12288 --use-singularity \
    --singularity-args '-B $PWD:/data --pwd /data' \
    --singularity-prefix=../00container
```

### Rule parameters

The `params` section of a rule is a good place for storing arguments
to commands used in the shell commands, tool versions, and other non-input
parameters. For example:

```python
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    threads: 4
    params: hisat = "-k 4"
    resources: mem_mb = 6144
    singularity:
        "library://wresch/classes/rnaseq:0.6"
    shell:
        """
        mkdir -p 02aln
        hisat2 {params.hisat} -x {input.idx} -U {input.fq} --threads {threads} \
          | samtools sort -T tmp/{wildcards.sample} -O BAM \
          > {output.bam}
        samtools index {output.bam}
        """
```

### Profiles

It can get pretty pretty cumbersome to re-type all the parameters for
snakemake. Snakemake profiles can set default values for command line
parameters.  A profile is a directory containing at a minimum a config.yaml
with keys corresponding to command line flags. A profile specified by
name on the command line is search for in `$HOME/.config/snakemake`. Alternatively
profiles can be specified by path. There is a simple profile available in
this directory which specifies all the required singularity settings so
we don't have to include them on the command line any more. Note that
profiles can be used to do much more sophisticated configuration.

```console
user@cn1234> cat myprofile/config.yaml
use-singularity: true
singularity-prefix: ../00container
singularity-args: '-B $PWD:/data --pwd /data'

user@cn1234> snakemake --cores 8 --profile ./myprofile
```
