## Parallelizing workflows across the cluster

**Goal:** Run the workflow in parallel by submitting some jobs to the cluster as batch jobs. As before,
use `Snakefile` as a starting point.

The changes required to the Snakefile in this case are minimal. We explicitly declare which rules
are to be executed on the same host as the main snakemake process even if submitting other jobs as
cluster batch jobs. At the top level of the Snakefile:

```python
localrules: all, clean
```

Snakemake can be taught to submit batch jobs by providing a template string that has access to
many of the properties of jobs (threads, resources, params) specified in the Snakefile, as well
as parameters defined in a cluster config file.  As for the general config file, the cluster
config file can be in yaml or json format. Here is an example that will work on biowulf:

```yaml
__default__:
    partition: quick
    time: 10
    extra: "--gres=lscratch:10"
hisat2:
    partition: norm
```

It makes sense to use this file for cluster specific settings such as partitions and local
scratch and runtimes. This file can contain a special `__default__` entry whose values will
be used for any rules that do not explicitly specify values for the keys in question. In this
example, hisat would be run on the norm partition with a walltime of 10 min and 10GB of
lscratch - the latter two from the `__default__` section since they are not defined for the
hisat2 rule.

When submitting to the cluster, a number of other of other options are important:

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

<div style="background-color:# FFFACD;">Please don't run snakemake workflows
on the login node, even if submitting jobs as batch jobs. Run the main process as a 
batch job itself or, if the workflow runs quickly enough, from an sinteractive
session</div>

So in our example (after adding the localrules declaration described earlier):

```console
user@cn1234> snakemake --jobs 6 --profile ./myprofile --cluster-config=cluster.yml \
    -w 120 -k --local-cores=10 --max-jobs-per-second 1 \
    --cluster="sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} \
        --time={cluster.time} {cluster.extra}"
```


And as before, most of these can be incorporated into a profile file for convenience:
```console
user@cn1234> cat myprofile/config.yaml
use-singularity: true
singularity-prefix: ../00container
singularity-args: '-B $PWD:/data,/lscratch/$SLURM_JOB_ID:/tmp --pwd /data'
max-jobs-per-second: 1
latency-wait: 120
keep-going: true
cluster: 'sbatch -c {threads} --mem={resources.mem_mb} --partition={cluster.partition} --time={cluster.time} {cluster.extra}'

user@cn1234> snakemake --profile ./myprofile --jobs 6 --cluster-config=cluster.yml --local-cores 10
```
