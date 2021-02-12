## Exercise 00 - Singularity refresher.

This is a brief reminder of how to use singularity containers. For more detail
see the materials for the companion
[Singularity class](https://github.com/NIH-HPC/Singularity-Tutorial) or the
[Singularity documentation](https://sylabs.io/guides/latest/user-guide/).

The container used in the exercises is `library://wresch/classes/rnaseq` and
should have also been downloaded locally if setup was successful. It can be
found in the `00container` directory in the root of the repository along with
its definition file:

```
REPO_ROOT/00container/
-rwxr-xr-x 1 user group 1.3G Feb 10 12:45 rnaseq.sif
-rw-r--r-- 1 user group 3.1K Feb 10 20:02 Singularity
```

It contains a number of tools that can be used in an RNA-Seq analysis (hisat2,
salmon, samtools, subread, R, ...). If you need to fetch the container locally you
can do so with

```console
user@cn1234> singularity pull library://wresch/classes/rnaseq
```

However, the setup script already took care of doing this, so there is not
need to run this command again.

### Recap of basic commands

If you are running this on the NIH HPC cluster, please make sure this and all
the following exercises are done within an interactive session with the
singularity and snakemake modules loaded:

```console
loginnode$ sinteractive --cpus-per-task=12 --gres=lscratch:20 --mem=24g
...
user@cn1234> module load singularity snakemake
[+] Loading singularity 3.7.1
[+] Loading snakemake 5.24.1
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

Singularity containers, in addition to being used with the singularity command,
are alse executable. When executed in this way, the runscript inside the
container is executed, which in our example prints a simple message.

```console
user@cn1234> # define a variable containing the full path to the container for connvenience
user@cn1234> container=$(cd .. && pwd)/00container/rnaseq.sif
user@cn1234> $container
------------------------------------------------------------
rnaseq - rnaseq pipeline tools version 0.5
------------------------------------------------------------
...
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

Note that this container does not have directories in the container
corresponding to all possible mount points on the NIH HPC systems. Therefore
singularity may not able to automatically bind mount the current working
directory inside the container. In that case it is also not able to set the
working directory inside the container to the same as outside. It is, however,
able to mount `/home` which then becomes the default working directory:

```console
user@cn1234> singularity shell $container
Singularity> pwd
/home/user
Singularity> exit
user@cn1234>
```

For the workflows in later examples to function, we will bind mount the local
directory under `/data` using `-B path_outside:path_inside` and change the
starting working directory to `/data` within the container with `--pwd` to
start in the same directory:

```console
user@cn1234> ls -lh
total 16K
-rw-r----- 1 user group 4.3K Feb  6 07:24 README.md
-rwxr-xr-x 1 user group  751 Feb  5 19:30 rnaseq
-rwxr-x--- 1 user group   95 Feb  5 14:13 test.sh

user@cn1234> singularity shell -B $PWD:/data   --pwd /data $container
#                 bind mount ---^----------^   ^---^-- starting directory
Singularity> pwd
/data
Singularity> ls -lh
total 16K
-rw-r----- 1 user group 4.3K Feb  6 07:24 README.md
-rwxr-xr-x 1 user group  751 Feb  5 19:30 rnaseq
-rwxr-x--- 1 user group   95 Feb  5 14:13 test.sh
```

### Wrapper script

`rnaseq` is an example script that uses the concepts introduced above to create a
wrapper for our Singularity container. It takes care to bind the current
directory under `/data` inside the container and passes all arguments on to a
bash shell inside the container via `singularity exec`. If lscratch exists, it
is mounted under `/tmp` inside the container. This wrapper can be used to
execute scripts or commands inside the container.

```console
user@cn1234> ./rnaseq --help
WARNING: Bind file source does not exist on host: /etc/resolv.conf
GNU bash, version 4.3.30(1)-release-(x86_64-pc-linux-gnu)
Usage:  bash [GNU long option] [option] ...
        bash [GNU long option] [option] script-file ...
...

user@cn1234> ./rnaseq -c "samtools --version"
samtools 1.10
```

Now, take for example the following `test.sh` script:

```bash
#! ./rnaseq
command -v hisat2
echo "-------------"
command -v salmon
echo "-------------"
```

If this script is executed by bash it will fail since none of these tools are
included in the path by default:

```console
user@cn1234> bash test.sh
--------------------
--------------------
```

However, execution with the wrapper works:

```console
user@cn1234> ./rnaseq test.sh
/opt/conda/envs/rnaseq/bin/hisat2
--------------------
/opt/salmon/1.4.0/bin/salmon
--------------------
```

And, more interestingly, if the script is made executable, it will be run from within
the container since the `#!` line used the wrapper script:

```console
user@cn1234> chmod +x test.sh
user@cn1234> ./test.sh
/opt/conda/envs/rnaseq/bin/hisat2
--------------------
/opt/salmon/1.4.0/bin/salmon
--------------------
```
