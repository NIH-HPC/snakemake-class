## Exercise 00

This is a short recap of how to use singularity containers. For more detail
see the materials for the companion
[Singularity class](https://github.com/NIH-HPC/Singularity-Tutorial) or the
[Singularity documentation](http://singularity.lbl.gov/quickstart).

The container used in the exercises and its definition file are located in the
`00container` directory in the root of the repository:

```
REPO_ROOT/00container/
|-- [wresch   1.0G]  NIH-HPC-snakemake-class-master-latest.simg
`-- [wresch   7.6K]  Singularity
```

It contains a number of tools that can be used in an RNA-Seq analysis (hisat2,
salmon, samtools, subread, R, ...). The repository is linked to the [Singularity
hub](https://www.singularity-hub.org/) which can generate new builds of containers
either automatically when changes are pushed to GitHub repo or manually on
request. The container can be pulled from Singularity hub with

```
singularity pull shub://NIH-HPC/snakemake-class
```

However, the setup script already took care of doing this, so there is not
need to run this command again.

#### Recap of basic commands

Make sure this and all the following exercises are done within an interactive
session with the singularity and snakemake modules loaded:

```console
loginnode$ sinteractive --cpus-per-task=12 --gres=lscratch:20
...
user@cn1234> module load singularity snakemake
[+] Loading singularity 2.4.1
[+] Loading snakemake 4.5.1
user@cn1234>
```

Singularity containers are executable. When executed in this way, the
runscript inside the container is executed, which in this case prints a simple
message.

```console
user@cn1234> # define a variable containing the full path to the container for connvenience
user@cn1234> container=$(cd .. && pwd)/00container/NIH-HPC-snakemake-class-master-latest.simg
user@cn1234> $container
------------------------------------------------------------
rnaseq - rnaseq pipeline tools version 0.3
------------------------------------------------------------
...
```

`singularity exec` is used to execute other programs inside the container

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
singularity is not able to automatically bind mount the current working
directory inside the container. It is able to mount `/home` which then becomes
the default working directory:

```console
user@cn1234> singularity shell $container
Singularity> pwd
/home/user
Singularity> exit
user@cn1234>
```

For the workflows in later examples to function, we will bind mount the local
directory under `/data` in the container and change the starting working
directory to `/data` within the container to start in the same directory:

```console
user@cn1234> ls -lh
total 16K
-rw-r----- 1 user group 4.3K Feb  6 07:24 README.md
-rwxr-xr-x 1 user group  751 Feb  5 19:30 rnaseq
-rwxr-x--- 1 user group   95 Feb  5 14:13 test.sh

user@cn1234> singularity shell -B $PWD:/data --pwd /data $container
Singularity> pwd
/data
Singularity> ls -lh
total 16K
-rw-r----- 1 user group 4.3K Feb  6 07:24 README.md
-rwxr-xr-x 1 user group  751 Feb  5 19:30 rnaseq
-rwxr-x--- 1 user group   95 Feb  5 14:13 test.sh
```

#### Wrapper script

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
samtools 1.5
```

Now, take for example the following `test.sh` script:

```bash
#! ./rnaseq
which hisat2
echo "-------------"
which salmon
echo "-------------"
```

If this script is executed by bash it will fail since none of these tools are
included in the path by default:

```console
user@cn1234> bash test.sh
which: no hisat2 in (...)
which: no salmon in (...)
```

However, execution with the wrapper works:

```console
user@cn1234> ./rnaseq test.sh
/opt/hisat2-2.1.0/hisat2
--------------------
/opt/Salmon-latest_linux_x86_64/bin/salmon
--------------------
```

And, more interestingly, if the script is made executable, it will be run from within
the container since the `#!` line used the wrapper script:

```console
user@cn1234> chmod +x test.sh
user@cn1234> ./test.sh
/opt/hisat2-2.1.0/hisat2
--------------------
/opt/Salmon-latest_linux_x86_64/bin/salmon
--------------------
```
