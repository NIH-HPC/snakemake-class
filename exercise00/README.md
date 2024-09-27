## Exercise 00 - Singularity refresher.

In this tutorial we will be using tools packaged into a singularity container
with our workflows. This helps portability and reproducibility of workflows.
Exercise 00 is a brief introduction to singularity. While not strictly required
for this class it is useful to understand the basics.

For more details see the official [Singularity documentation](https://sylabs.io/guides/latest/user-guide/)
and our [Singularity overview](https://hpc.nih.gov/apps/singularity.html).

The container used in the exercises is `library://wresch/classes/rnaseq` and
should have locally cached during setup. After runnign the setup script it
can be found in the `00container` directory in the root of the repository along with
its definition file:

```
REPO_ROOT/00container/
-rwxr-xr-x 1 user group 1.3G Feb 10 12:45 2354d2ff28bcf0b42c57fae398b4c9b5.simg
lrwxrwxrwx 1 user group   37 Sep 26 14:31 rnaseq.sif -> 2354d2ff28bcf0b42c57fae398b4c9b5.simg
-rw-r--r-- 1 user group 3.1K Feb 10 20:02 rnaseq.def
```

It contains a number of tools that can be used in an RNA-Seq analysis (hisat2,
salmon, samtools, subread, R, ...). If you wanted to manually fetch this
container you could do so with

```console
user@cn1234> singularity pull library://wresch/classes/rnaseq
```

### What is singularity

> SingularityCE is a container platform. It allows you to create and run
> containers that package up pieces of software in a way that is portable and
> reproducible. You can build a container using SingularityCE on your laptop, and
> then run it on many of the largest HPC clusters in the world, local university
> or company clusters, a single server, in the cloud, or on a workstation down
> the hall. Your container is a single file, and you don’t have to worry about
> how to install all the software you need on each different operating system.
> 
> _From_ [SingularityCE introduction](https://docs.sylabs.io/guides/latest/user-guide/introduction.html)

Containers in general package software an all its dependencies into a package
that can be executed in a more or less isolated environment. Singularity and
the closely related Apptainer in particular are an implementation of container
technology geared towards scientific computing.

### Recap of basic singularity commands

If you are running this on the NIH HPC cluster, please make sure this and all
the following exercises are done within an interactive session with the
singularity and snakemake modules loaded:

```console
loginnode$ sinteractive --cpus-per-task=12 --gres=lscratch:20 --mem=24g
...
user@cn1234> module load singularity snakemake
[+] Loading singularity  4.0.3  on cn4301
[+] Loading snakemake  7.32.4

user@cn1234>
```

singularity allows you to create and use single file containers. You can run
programs inside a container transparently: IO can be redirected; arguments
passed, and files accessed; Ports on the host system from within a container
and a user retains their identity inside the container.

The singularity command includes many subcommands and options. use

```console
user@cn1234> singularity help
```

to get more information about this tool which is the main way of interacting
with singularity containers.

Singularity container files, in addition to being used with the singularity command,
can also be made executable. When executed in this way, the runscript inside the
container is executed, which in our example prints a simple message.

```console
user@cn1234> # define a variable containing the full path to the container for connvenience
user@cn1234> container=$(cd .. && pwd)/00container/rnaseq.sif
user@cn1234> $container
------------------------------------------------------------
rnaseq - rnaseq pipeline tools version 0.8
------------------------------------------------------------

This container encapsulates tools for RNA-Seq analysis.
It is intended for creating reproducible pipelines.

```

`singularity exec` is used to execute programs inside the container

```console
user@cn1234> singularity exec $container samtools index
Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
Options:
  -b       Generate BAI-format index for BAM files [default]
  -c       Generate CSI-format index for BAM files
  -m INT   Set minimum interval size for CSI indices to 2^INT [14]
  -@ INT   Sets the number of threads [none]
```

By default, programs running inside a container have minimal access
to the file systems outside the container. Let's start a shell inside our
example container and try to access the /fdb directory, for example:

```console
user@cn1234> singularity shell $container
Singularity> ls -lh /fdb                           
ls: cannot access '/fdb': No such file or directory
```

The only external paths accessible are your home directory and the current
directory (`/data/username/snakemake-class/exercise00` in this example):

```console
Singularity> pwd
/data/username/snakemake-class/exercise00
Singularity> cd ~
Singularity> pwd
/home/username
Singularity> ls
Desktop bin temp ...
Singularity> exit
user@cn1234>
```

However, other paths can be made visible to processes inside the container via
what is called a bind mount. Bind mounts can be define on the command line with `-B`
or by setting the `SINGULARITY_BINDPATH` variable. For example

```console
user@cn1234> singularity shell -B /fdb/STAR_indices $container
Singularity> ls -l /fdb
total 0                                                 
drwxr-xr-x 2 user group 4096 Apr 12 11:35 STAR_indices
Singularity> exit

user@cn1234> ## use outside_path:inside_path to make paths visible under a different path
                to containerized processes
user@cn1234> singularity shell -B /fdb/STAR_indices:/star,/fdb/salmon:/salmon $container
Singularity> ls -l /fdb
ls: cannot access '/fdb': No such file or directory
Singularity> ls -l /star
total 58
...
drwxrwxr-x 2 user group  4096 Sep 26  2018 2.5.4a
...
```

We provide a convenient script you can source to set up a bindpath
variable that will make all our file systems accessible at the same paths
inside the container

```console
user@cn1234> source /usr/local/current/singularity/app_conf/sing_binds
user@cn1234> echo $SINGULARITY_BINDPATH
/gs10,/gs11,/gs12,/vf,/spin1,/data,/fdb,/gpfs,/lscratch
user@cn1234> singularity shell $container
Singularity> ls /fdb
00_TO_BE_DELETED   T2T       dbscSNV       humann     purge_haplotigs
...
Singularity> exit
```

For the workflows in this tutorial we will use the `SINGULARITY_BINDPATH`
defined by our helper script so that paths inside and outside the container
will be identical. If you only use relative paths this is not in fact
necessary.

