## Exercise 00

This is a short recap of how to use singularity containers. For more detail
see the materials for the companion
[Singularity class](https://github.com/NIH-HPC/Singularity-Tutorial) or the
[Singularity documentation](http://singularity.lbl.gov/quickstart).

The container used in the exercises and its definition file are located in the
`00container` directory in the root of the repository:

```
$REPO_ROOT/00container/
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

Make sure this and all the following exercises are done within an interactive
session with the singularity and snakemake modules loaded:

```
loginnode$ sinteractive --cpus-per-task=12 --gres=lscratch:20
...
cn1234$ module load singularity snakemake
[+] Loading singularity 2.4.1
[+] Loading snakemake 4.5.1
cn1234$
```

Singularity containers are executable. When executed in this way, the
runscript inside the container is executed, which in this case prints a simple
message.

```
cn1234$ # define a variable containing the full path to the container for connvenience
cn1234$ container=$(cd .. && pwd)/00container/NIH-HPC-snakemake-class-master-latest.simg
cn1234$ $container
------------------------------------------------------------
rnaseq - rnaseq pipeline tools version 0.3
------------------------------------------------------------
...
```

`singularity exec` is used to execute other programs inside the container
```
cn1234$ singularity exec $container samtools index
Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
Options:
  -b       Generate BAI-format index for BAM files [default]
  -c       Generate CSI-format index for BAM files
  -m INT   Set minimum interval size for CSI indices to 2^INT [14]
  -@ INT   Sets the number of threads [none]
```

`rnaseq` is a wrapper script present in each exercise directory. It takes care
to bind the current directory under `/data` inside the container and passes all
arguments on to a bash shell inside the container via `singularity exec`. That
means it can be used to execute scripts or commands inside the container.

```
cn1234$ ./rnaseq --help
WARNING: Bind file source does not exist on host: /etc/resolv.conf
GNU bash, version 4.3.30(1)-release-(x86_64-pc-linux-gnu)
Usage:  bash [GNU long option] [option] ...
        bash [GNU long option] [option] script-file ...
...

cn1234$ ./rnaseq -c "samtools --version"
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

```
cn1234$ bash test.sh
which: no hisat2 in (...)
which: no salmon in (...)
```

However, execution with the wrapper works:
```
cn1234$ ./rnaseq test.sh
/opt/hisat2-2.1.0/hisat2
--------------------
/opt/Salmon-latest_linux_x86_64/bin/salmon
--------------------
```

And, more interestingly, if the script is made executable, it will be run from within
the container since the `#!` line used the wrapper script:

```
cn1234$ chmod +x test.sh
cn1234$ ./test.sh
/opt/hisat2-2.1.0/hisat2
--------------------
/opt/Salmon-latest_linux_x86_64/bin/salmon
--------------------
```
