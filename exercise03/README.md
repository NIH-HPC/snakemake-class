## Add a new rule to count RNA-Seq reads per gene

**Goal:** Add a new rule that summarizes the count of RNA-Seq reads per gene
for all 6 samples. The starting point, as before, is `Snakefile` and an example
of a finished workflow is provided in `Snakefile.finished`


### Add the featureCounts rule

`featureCounts` from the subread package is used to count the number of reads
from each alignment mapping to each gene. It's **inputs** are a bam file and
the annotation file `00ref/R64-1-1.genes.gtf`. It's **output** is a per sample
count file in the `04count` directory. It's expected to use 4GB of memory and
run 2 threads.

Here is a possible implementation of this rule:

```python
rule count:
    input: bam = "02aln/{sample}.bam",
           annot = "00ref/R64-1-1.genes.gtf"
    output: "04count/{sample}"
    threads: 2
    resources: mem_mb = 4096
    singularity:
        "library://wresch/classes/rnaseq:0.5"
    shell:
        """
        featureCounts -a {input.annot} -o {output} \
                -T {threads} --minOverlap 10 {input.bam}
        """
```

Now, since the final output of the pipeline as it is are the count files, the
`all` rule has to be changed to request the generation of the count files.
Alignment files don't have to be specified any more b/c snakemake will
automatically determine that alignments are required to genrate count files.
The new `all` rule would therefore be

```python
rule all:
    input: "04count/ERR458495",
           "04count/ERR458502"
```

### Extend the pipeline to include all 6 samples

Snakefiles are essentially python and arbitrary python code can be used in
many places. Let's take advantage of this and extend the workflow to include
all fastq files present in the `00fastq` directory:

```python
import os.path
from glob import glob

# use the glob function to get all fastq files from the 00fastq directory
# extract the sample name from each path (assumes that there is one fastq per sample)
samples = [os.path.basename(a).replace(".fastq.gz", "") for a in glob("00fastq/*.fastq.gz")]

# the all rule creates a list of all count files that sould be generated based on
# the list of samples
rule all:
    input: expand("04count/{s}", s=samples)
```

Now the workflow can be run on all samples in parallel. Again using a profile that
sets the correct options for using the singularity container.

```console
user@cn1234> snakemake --profile ./myprofile --cores 12 --resources mem_mb=12288
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 12
Rules claiming more threads will be scaled down.
Provided resources: mem_mb=12288
Job counts:
        count   jobs
        1       all
        6       count
        6       hisat2
        13
...

