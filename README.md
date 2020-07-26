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
./
├── 0_Preparation
│   ├── pombe_ASM294v1_18_toplevel.dict
│   ├── pombe_ASM294v1_18_toplevel.fa
│   ├── pombe_ASM294v1_18_toplevel.fa.amb
│   ├── pombe_ASM294v1_18_toplevel.fa.ann
│   ├── pombe_ASM294v1_18_toplevel.fa.bwt
│   ├── pombe_ASM294v1_18_toplevel.fa.fai
│   ├── pombe_ASM294v1_18_toplevel.fa.pac
│   ├── pombe_ASM294v1_18_toplevel.fa.sa
│   └── Schizosaccharomyces_pombe.ASM294v1.18_noSPBC3F6.03.gtf
├── 1_Files
│   ├── 0_RawData
│   │   ├── NSK-0912-176_AK9183-N505_1.fastq
│   │   └── NSK-0912-176_AK9183-N505_2.fastq
│   ├── 1_Fastp
│   │   ├── NSK-0912-176_AK9183-N505_1.fastp.fastq
│   │   ├── NSK-0912-176_AK9183-N505_2.fastp.fastq
│   │   ├── NSK-0912-176_AK9183-N505.html
│   │   └── NSK-0912-176_AK9183-N505.json
│   ├── 2_BWA
│   │   ├── NSK-0912-176_AK9183-N505
│   │   └── NSK-0912-176_AK9183-N505.bwa.sam
│   ├── 3_ProcessedBAM
│   │   ├── NSK-0912-176_AK9183-N505.bwa.srt.bam
│   │   ├── NSK-0912-176_AK9183-N505.bwa.srt.mark_dup.bai
│   │   ├── NSK-0912-176_AK9183-N505.bwa.srt.mark_dup.bam
│   │   └── NSK-0912-176_AK9183-N505.marked_dup_metrics.txt
│   └── 4_Variants_Samtools
│       └── NSK-0912-176_AK9183-N505.Samtools.vcf
├── 2_Summary
└── scripts
```