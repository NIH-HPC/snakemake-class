
Setup - please complete before class
================================================================================


After cloning this repository, a number of steps (caching the singularity
container and data, creating reference files, ...) have to be carried out
before starting class. This requires `singularity`, `snakemake` and `git`
and has to be executed on a compute node.



On the NIH HPC systems, start an interactive session, load the snakemake and
singularity modules, and clone this repository:

```console
user@headnode> sinteractive --cpus-per-task=12 --mem=24g --gres=lscratch:20
...
user@cn1234> ## change to a suitable directory somewhere in /data
user@cn1234> cd /data/$USER
user@cn1234> module load git
user@cn1234> git clone https://github.com/NIH-HPC/snakemake-class.git
user@cn1234> cd snakemake-class
```

Setup will be different on other systems.

When you are ready run the setup script to fetch data and create all files
necessary for the exercies.

```
user@cn1234> ./setup.sh
...
+------------------------------------------------------------------------------+
|                                                                              |
|                Class materials have been set up successfully                 |
|                                                                              |
+------------------------------------------------------------------------------+
```

The `setup.smk` in the root directory of the repo takes care of all the setup
required for class. It also serves as another example for a Snakemake workflow.

## Exercises

* [Exercise 0 - Singularity refresher](/exercise00/)
* [Exercise 1 - Basic snakemake rules](/exercise01/)
* [Exercise 2 - Parallelizing pipelines](/exercise02/)
* [Exercise 3 - Adding more rules](/exercise03/)
* [Exercise 4 - Running on a HPC cluster](/exercise04/)
* [Exercise 5 - Visualizing workflows](/exercise05/)
