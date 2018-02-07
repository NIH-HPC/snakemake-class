[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/558) 

Setup - please complete before class
================================================================================


After cloning this repository, a number of steps (downloading the singularity
container and data, ...) have to be carried out to complete before starting
class. This requires that `singularity` and `snakemake` are on the path
for the setup. Here is an overview of the complete setup:

On the NIH HPC systems start an interactive session and clone the repository:
```
user@headnode> sinteractive --cpus-per-task=12 --mem=24g --gres=lscratch:20
...
user@cn1234> module load singularity snakemake hisat2
user@cn1234> cd /data/$USER # or whereever you'd like the class directory to be
user@cn1234> git clone https://github.com/NIH-HPC/snakemake-class.git
user@cn1234> cd snakemake-class
```

Then run the setup
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
