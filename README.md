[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/558) 

Setup
================================================================================

After cloning this repository, a number of steps (downloading the singularity
container and data, ...) have to be carried out to complete before starting
class. This requires `singularity` and `snakemake`. Here is an overview of the
complete setup:

On the NIH HPC systems start an interactive session and clone the repository:
```bash
$ sinteractive --cpus-per-task=10 --gres=lscratch:20
...
$ module load singularity snakemake
$ cd /data/$USER
$ git clone https://github.com/NIH-HPC/snakemake-class.git
$ cd snakemake-class
```

Then run the setup
```bash
$ snakemake setup
...
```

