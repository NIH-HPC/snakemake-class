## Parallelizing workflows across the cluster

**Goal:** Run the workflow in parallel by submitting some jobs to the cluster
as batch jobs. As before, use `Snakefile` as a starting point.

The changes required to the Snakefile in this case are minimal. We explicitly
declare which rules are to be executed on the same host as the main snakemake
process even if submitting other jobs as cluster batch jobs. At the top level
of the Snakefile:

```python
localrules: all, clean
```

Snakemake can be taught to submit batch jobs by providing a template string
that has access to many of the properties of jobs (threads, resources, params)
specified in the Snakefile, as well as parameters defined in a cluster config
file. This can be done with plain snakemake. However, there is a Biowulf-specific
[profile](https://github.com/NIH-HPC/snakemake_profile) that takes care of of job submission
and job status checking without taxing the slurm scheduler. Rather than refining
the profile from the previous exercise, we will switch to the Biowulf profile.
Instead of specifying the profile manually every time we run snakemake
we'll set an environment variable

```console
user@cn1234> export SNAKEMAKE_PROFILE="$(cd .. && pwd)/bwprofile"
user@cn1234> cat $SNAKEMAKE_PROFILE/config.yaml
restart-times: 0
jobscript: slurm_jobscript.sh
cluster: bw_submit.py
cluster-status: bw_status.py
cluster-cancel: scancel
max-jobs-per-second: 1
max-status-checks-per-second: 1
local-cores: 4
latency-wait: 240
jobs: 50
keep-going: True
rerun-incomplete: True
# note that the space before the -- is necessary - otherwise cluster
# execution failes as --cleanenv gets interpreted as a new option since
# snakemake passes this on as --singularity-args "--cleanenv" - i.e. without
# a `=`
singularity-args: " --cleanenv"
use-singularity: true
singularity-prefix: /data/user/snakemake-class/00container
```

The last line will be specific to your class directory. Use `snakemake --help`
to understand the different options specified in the profile.
In brief, here are some of the most relevant options:

  - `-k`, `--keep-going`: By default, snakemake will quit if a job fails (after waiting
    for running jobs to finish. `-k` will make snakemake continue with independent jobs.
  - `-w`, `--latency-wait`, `--output-wait`: The amount of time snakemake will wait for
    output files to appear after a job has finished. This defaults to a low 5s. On the
    shared file systems latency may be higher. Raising it to 120s is a bit excessive, but
    it doesn't really hurt too much.
  - `--local-cores`: The number of CPUs available for local rules
  - `--max-jobs-per-second`: Max numbers of jobs to submit per second. Please be kind
    to the batch scheduler.
  - `--cluster`: The template string used to submit each (non local) job.
  - `--jobs`: The number of jobs to run concurrently.
  - `--cluster-config`: The cluster config file

:information_source: Please **do not run snakemake workflows on the login node**, 
even if submitting jobs as batch jobs. Run the main process as a 
batch job itself or, if the workflow runs quickly enough, from an sinteractive
session

So in our example (after adding the `localrules` declaration described earlier and
setting `$SNAKEMAKE_PROFILE`):

```console
user@cn1234> snakemake 
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cluster nodes: 50
Job stats:
job       count
------  -------
all           1
count         6
hisat2        6
total        13

Select jobs to execute...

[Thu Sep 26 18:48:20 2024]
rule hisat2:
    input: 00fastq/ERR458495.fastq.gz, 00ref/hisat_index/R64-1-1
    output: 02aln/ERR458495.bam, 02aln/ERR458495.bam.bai
    jobid: 8
    reason: Missing output files: 02aln/ERR458495.bam
    wildcards: sample=ERR458495
    threads: 4
    resources: mem_mb=6144, mem_mib=5860, disk_mb=1000, disk_mib=954, tmpdir=<TBD>

hisat2: submission command "sbatch --cpus-per-task=4 --mem=6144 --time=120 --gres=lscratch:1 --output=logs/hisat2-%j.out --partition=quick /spin1/users/wresch/code/class_materials/snakemake-class/exercise05/.snakemake/tmp.cxnvsixj/snakejob.hisat2.8.sh
Submitted job 8 with external jobid '36450266'.

...
```


