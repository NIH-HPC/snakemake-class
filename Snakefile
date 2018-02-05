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
           "00ref/R64-1-1.fa", 
           "00ref/hisat_index/R64-1-1", 
           "00ref/R64-1-1.cdna_nc.fa", 
           "00ref/R64-1-1.genes.gtf",
           "00ref/ref.yml",
           "00ref/R64-1-1.tran2gene.tsv",

###
### exercises and container
###
rule fetch_container:
    output: "00container/rnaseq.simg"
    shell: 
        """
        module load singularity
        # for some reason pulling to subdir does not work for other users...
        cd 00container
        singularity pull -n rnaseq.simg shub://NIH-HPC/snakemake-class
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
rule mkdir_00ref:
    output: "00ref"
    shell: "mkdir 00ref"

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
    input: "00fastq/samples.yml", "00ref/ref.yml"
    output: "00setup/config.yml"
    shell:
        """
        echo "container: 00config/rnaseq.img" > {output}
        echo "wrapper: rnaseq" >> {output}
        cat {input} >> {output}
        """


###
### reference data
###

ENSEMBL_RELEASE = 91
ENSEMBL_URL = "ftp://ftp.ensembl.org/pub/release-{}".format(ENSEMBL_RELEASE)

genome_build = "R64-1-1"


rule fetch_genome:
    output: "00ref/R64-1-1.fa"
    shell:
        """
        wget {ENSEMBL_URL}/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
        gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
        mv Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa {output}
        """

rule fetch_transcriptome:
    output: "00ref/R64-1-1.cdna_nc.fa"
    shell:
        """
        wget {ENSEMBL_URL}/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
        wget {ENSEMBL_URL}/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz
        gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
        gunzip Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz
        cat Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa \
            Saccharomyces_cerevisiae.R64-1-1.ncrna.fa \
              | awk '/^>/ {{NF=1}} {{print}}' > {output}
        rm Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa \
           Saccharomyces_cerevisiae.R64-1-1.ncrna.fa
        """

rule fetch_gtf:
    output: "00ref/R64-1-1.genes.gtf"
    shell:
        """
        wget "{ENSEMBL_URL}/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.{ENSEMBL_RELEASE}.gtf.gz"
        gunzip Saccharomyces_cerevisiae.R64-1-1.{ENSEMBL_RELEASE}.gtf.gz
        mv Saccharomyces_cerevisiae.R64-1-1.{ENSEMBL_RELEASE}.gtf {output}
        """

rule gtf2bed12:
    input: "00ref/R64-1-1.genes.gtf"
    output: "00ref/R64-1-1.genes.bed12"
    run:
        transcripts = {}
        i = 0
        for line in open(input[0]):
            if line.startswith('#'):
                continue
            chrom, _, feat, s, e, score, strand, _, attr = line.split('\t')
            if feat not in ('transcript', 'exon'):
                continue
            s = int(s) - 1
            e = int(e)
            assert s < e
            tid = None
            for a in attr.split(';'):
                if a.strip().startswith('transcript_id'):
                    tid = a.split('"')[1]
                    break
            if tid is None:
                raise ValueError
            if feat == 'transcript':
                i += 1
                transcripts[tid] = {'c': chrom, 's':s, 'e':e, 
                        'strand': strand, 'exons':[], 'tid': tid, 'i': i}
            else:
                transcripts[tid]['exons'].append((s, e))
        with open(output[0], "w") as of:
            for tid, t in sorted(transcripts.items(), key = lambda x: x[1]['i']):
                of.write("{c}\t{s}\t{e}\t{tid}\t0\t{strand}".format(**t))
                of.write("\t{s}\t{e}\t0,0,0\t{en}".format(en=len(t['exons']), **t))
                exons = sorted(t['exons'])
                of.write("\t{bsz}\t{bs}\n".format(bsz=",".join(str(e-s) for s,e in exons),
                    bs=",".join(str(s) for s,_ in exons)))

rule make_transcript_gene_map:
    input: "00ref/R64-1-1.genes.gtf"
    output: "00ref/R64-1-1.tran2gene.tsv"
    threads: 1
    run:
        transcripts = []
        for line in open(input[0]):
            if line.startswith('#'):
                continue
            chrom, _, feat, s, e, score, strand, _, attr = line.split('\t')
            if feat != 'transcript':
                continue
            tid = None
            gid = None
            for a in attr.split(';'):
                if a.strip().startswith('transcript_id'):
                    tid = a.split('"')[1]
                if a.strip().startswith('gene_id'):
                    gid = a.split('"')[1]
            if tid is None or gid is None:
                raise ValueError
            transcripts.append((tid, gid))
        with open(output[0], "w") as of:
            for tid, gid in transcripts:
                of.write("{}\t{}\n".format(tid, gid))


rule ref_yaml:
    output: "00ref/ref.yml"
    shell:
        """
        echo "reference:" > {output}
        echo "  ensembl_ver: {ENSEMBL_RELEASE}" >> {output}
        echo "  genome_build: R64-1-1" >> {output}
        echo "  genome_file: R64-1-1.fa" >> {output}
        echo "  cdna_file: R64-1-1.cdna.fa" >> {output}
        """


rule make_hisat_index:
    input:  "00ref/R64-1-1.fa"
    output: idxf1 = "00ref/hisat_index/R64-1-1.1.ht2", 
            name = "00ref/hisat_index/R64-1-1"
    threads: 2
    shell:
        """
        ../rnaseq hisat2-build -p {threads} {input} {output.name} \
                && touch {output.name}
        """
