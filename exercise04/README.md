## Configuration of workflows

**Goal:** Extract sample information and the path to the hisat index from
a configuration file in yaml. As in the other exercises, the starting
file is `Snakemake` and the final product is `Snakemake.finished`.

Snakemake workflows can make use of configuration files in [yaml](yaml.org)
or [json](https://www.json.org/) format. The configuration file is specified
at the top level of the Snakefile. For example:

```python
configfile: "config.yml"

rule all:
    input: expand("04count/{s}", s=samples)
```

The config file is parsed and made available as a global dictionary named
`config`. So, given the following config file

```yaml
samples:
  ERR458502:
    gt: SNF2
    rep: 1
  ERR458509:
    gt: SNF2
    rep: 2
  ERR458516:
    gt: SNF2
    rep: 3
  ERR458495:
    gt: WT
    rep: 1
  ERR458880:
    gt: WT
    rep: 2
  ERR458887:
    gt: WT
    rep: 3
reference:
  ensembl_ver: 88
  genome_build: R64-1-1
  hisat_index: 00ref/hisat_index/R64-1-1
  genome_file: 00ref/R64-1-1.fa
  cdna_file: 00ref/R64-1-1.cdna_nc.fa
```

the `config` dict would look like this:

```python
{'reference': {'cdna_file': '00ref/R64-1-1.cdna_nc.fa',
               'ensembl_ver': 88,
               'genome_build': 'R64-1-1',
               'genome_file': '00ref/R64-1-1.fa',
               'hisat_index': '00ref/hisat_index/R64-1-1'},
 'samples': {'ERR458495': {'gt': 'WT', 'rep': 1},
             'ERR458502': {'gt': 'SNF2', 'rep': 1},
             'ERR458509': {'gt': 'SNF2', 'rep': 2},
             'ERR458516': {'gt': 'SNF2', 'rep': 3},
             'ERR458880': {'gt': 'WT', 'rep': 2},
             'ERR458887': {'gt': 'WT', 'rep': 3}}}
```

The sample list can be extracted like so:
```python
configfile: "config.yml"
samples = config["samples"].keys()

rule all:
    input: expand("04count/{s}", s=samples)
```

and the hisat index in the hisat rule like this:
```python
rule hisat2:                                       
    input: fq = "00fastq/{sample}.fastq.gz",       
           idx = config["reference"]["hisat_index"]
    ...
```
