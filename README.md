[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/558) 

Setup - please complete before class
================================================================================


After cloning this repository, a number of steps (downloading the singularity
container and data, ...) have to be carried out to complete before starting
class. This requires that `singularity` and `snakemake` are on the path
for the setup. Here is an overview of the complete setup:

On the NIH HPC systems start an interactive session and clone the repository:
```
loginnode$ sinteractive --cpus-per-task=10 --gres=lscratch:20
...
computenode$ module load singularity snakemake hisat2
computenode$ cd /data/$USER # or whereever you'd like the class directory to be
computenode$ git clone https://github.com/NIH-HPC/snakemake-class.git
computenode$ cd snakemake-class
```

Then run the setup
```
computenode$ snakemake setup
computenode$ snakemake --use-singularity \
    --singularity-args '-B $PWD:/data --pwd /data' \
    --singularity-prefix=00container setup
...
```

The `Snakefile` in the root directory of the repo takes care of all the setup
required for class. It also serves as another example for a Snakemake workflow.
