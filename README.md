# Variant Calling
This project is to use different software for variant calling  

## Introduction
In the whole genome sequencing analysis, variant calling is the process that finds the mutations on the genome by comparing the sequence information between the reference genome and the sequence reads. The abilities of variant calling software are different, therefore choosing the proper software is important during the analysis.

Here we tried samtools/bcftools, GATK, Deepvaraint to analyze the sequencing data.


## Variant Callers
| Caller | WebSite |
| :- | :- |
| Samtools/bcftools | http://www.htslib.org/ |
| GATK | https://gatk.broadinstitute.org/hc/en-us |
| Deepvariant | https://github.com/google/deepvariant |


# Usage
```shell
Usage: python WGS_Analysis_Pipeline.py [options]

    options:
        -h/--help
        -P/--project path to the project, such as /data/yangyusheng/projectName
        -R/--reference path to the reference genome
        -A/--annotation path to the annotation
        -D/--data raw data, folder stores the raw data, like /data/yangyusheng/
        -S/--suffix raw data suffix, the suffix of *.fastq is .fastq
        -C/--caller Samtools, GATK, Deepvariant
        -B/--BQSR
```
This script will create the folders with the following structure.
```

```