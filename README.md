
Setup - please complete before class
================================================================================


After cloning this repository, a number of steps (downloading the singularity
container and data, ...) have to be carried out before starting class. This
requires that `singularity`, `snakemake` and `git`.

On the NIH HPC systems, start an interactive session, load the snakemake and
singularity modules, and clone this repository:

```
user@headnode> sinteractive --cpus-per-task=12 --mem=24g --gres=lscratch:20
...
user@cn1234> cd /data/$USER # or whereever you'd like the class directory to be
user@cn1234> git clone https://github.com/NIH-HPC/snakemake-class.git
user@cn1234> module load snakemake singularity
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

The `Snakefile` in the root directory of the repo takes care of all the setup
required for class. It also serves as another example for a Snakemake workflow.

## Exercises

* [Exercise 0 - Singularity refresher](/exercise00/README.md)
* [Exercise 1 - Basic snakemake rules](/exercise01/README.md)
* [Exercise 2 - Parallelizing pipelines](/exercise02/README.md)
* [Exercise 3 - Adding more rules](/exercise03/README.md)
* [Exercise 4 - Running on a HPC cluster](/exercise04/README.md)
* [Exercise 5 - Visualizing workflows](/exercise05/README.md)
