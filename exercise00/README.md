## Exercise 00

This is a short recap of how to use singularity containers. The `00container`
directory contains a singularity container with a number of tools for RNAseq and
the definition file used to build the container. The repository is linked to
singularity hub, so the container gets built automatically and can be pulled with 

```
singularity pull shub://NIH-HPC/snakemake-class
```

`rnaseq` is a wrapper script for the container that bind mounts the current
working dir as /data inside the container. It can be used to execute commands
or scripts inside the container.

Load the singularity module. Note that this is only possible on a compute node
i.e. in a batch job or an interactive session:

```
$ module load singularity
[+] Loading singularity 2.4.1 on cn3140
```

Singularity containers are executable. When executed in this way, the
runscript inside the container is executed:
```
$ 00container/rnaseq.simg
```

`singularity exec` is used to execute other programs inside the container
```
$ singularity exec 00container/rnaseq.simg samtools index
```

The `rnaseq` wrapper script allows execution of commands or scripts
from withing the container and takes care to set up bind paths.
```
$ ./rnaseq -c "samtools"
```

Take for example the following `test.sh` script:
```bash
#! ./rnaseq
which hisat2
echo "-------------"
which salmon
echo "-------------"
```

It uses the wrapper script in it's shebang line. If this script is executed by
bash it will fail since none of these tools are included in the path by default:

```
$ bash test.sh
which: no hisat2 in (...)
which: no salmon in (...)
```

However, execution with the wrapper works:
```
$ ./rnaseq test.sh
/opt/hisat2-2.1.0/hisat2
--------------------
/opt/Salmon-latest_linux_x86_64/bin/salmon
--------------------
```

And, more interestingly, if the script is made executable, it will be run from within
the container since the `#!` line used the wrapper script:
```
$ chmod +x test.sh
$ ./test.sh
/opt/hisat2-2.1.0/hisat2
--------------------
/opt/Salmon-latest_linux_x86_64/bin/salmon
--------------------
```

Finally, the container can also be run directly from Singularity hub:
```
$ singularity exec shub://NIH-HPC/snakemake-class salmon --version
```
