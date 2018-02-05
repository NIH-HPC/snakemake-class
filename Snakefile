# vim: set ft=python:

rule setup:
    input: "00container/rnaseq.simg",
           expand("exercise{n:02d}/00ref", n=range(1)),
           expand("exercise{n:02d}/00fastq", n=range(1)),
           expand("exercise{n:02d}/rnaseq", n=range(1)),
           "exercise00/00container"

rule fetch_container:
    output: "00container/rnaseq.simg"
    shell: 
        """
        module load singularity
        singularity pull -n {output} shub://NIH-HPC/snakemake-class
        """

rule link_dirs:
    input: "00{dir}"
    output: "{ex}/00{dir}"
    shell: "cd {wildcards.ex} && ln -s ../{input}"

rule wrapper:
    input: "00setup/rnaseq"
    output: "{ex}/rnaseq"
    shell: "cp {input} {output}"
