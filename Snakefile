# vim: set ft=python:

SAMPLES = ["ERR458495",
        "ERR458502",
        "ERR458509",
        "ERR458516",
        "ERR458880",
        "ERR458887"]

rule setup:
    input: "00container/rnaseq.simg",
           "00fastq/ERP004763_sample_table.tsv",
           "00fastq/samples.yml",
           "00setup/config.yml",
           expand("00fastq/{sample}.fastq.gz", sample=SAMPLES),
           expand("exercise{n:02d}/00ref", n=range(1)),
           expand("exercise{n:02d}/00fastq", n=range(1)),
           expand("exercise{n:02d}/rnaseq", n=range(1)),
           "exercise00/00container",
rule all:
    input: "samples.yml",
           expand("{sample}.fastq.gz", sample=SAMPLES)

###
### exercises and container
###
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


###
### data
###

rule mkdir_00fastq:
    output: "00fastq"
    shell: "mkdir 00fastq"

rule fetch_sample_desc:
    """fetch the sample description table; use local repo if possible"""
    output: "00fastq/ERP004763_sample_table.tsv"
    shell:
        """
        if [[ -f /data/classes/snakemake/{output} ]]; then
            cp /data/classes/snakemake/{output} {output}
        else
            wget -O {output} \
              https://ndownloader.figshare.com/files/2194841
        fi
        """

rule sample_table:
    input: "00fastq/ERP004763_sample_table.tsv"
    output: "00fastq/samples.yml"
    shell:
        """
        echo "samples:" > {output}
        pattern="$(echo {SAMPLES} | tr -s ' ' '|')"
        egrep "$pattern" {input} \
                | sort -k3,4 \
                | awk '{{printf("  %s:\\n    gt: %s\\n    rep: %s\\n", $1, $3, $4)}}' \
                >> {output}
        """

rule fetch_fastq:
    """fetch fastq; use local repo if possible"""
    output: "00fastq/{sample}.fastq.gz"
    shell:
        """
        if [[ -f /data/classes/snakemake/{output} ]]; then
            cp /data/classes/snakemake/{output} {output}
        else
            url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
            sample="{wildcards.sample}"
            wget -O {output} "$url/${{sample:0:6}}/$sample/${{sample}}.fastq.gz"
        fi
        """

rule config_yml:
    input: "00fastq/samples.yml"
    output: "00setup/config.yml"
    run:
        with open(output[0], "w") as ofh:
            ofh.write("""container: 00config/rnaseq.img
wrapper: rnaseq
{samples}
reference:
  ensembl_ver: 88
  genome_build: R64-1-1
  hisat_index: 00ref/hisat_index/R64-1-1
  genome_file: 00ref/R64-1-1.fa
  cdna_file: 00ref/R64-1-1.cdna_nc.fa
""".format(samples=open(input[0]).read().strip()))
