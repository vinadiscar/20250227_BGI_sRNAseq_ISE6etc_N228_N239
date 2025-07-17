# 20250227_BGI_sRNAseq_ISE6etc_N228_N239
- This repository aims to analyze small RNA seq data. Libraries are in N228 -N239 from project disc /Volumes/okamura-lab/raw_sequencing_data/20250227_BGI_sRNAseq_ISE6etc

## I. Preprocess
### 1a. fastp
```sh
fastp.md
fastp_SE_BGI_sRNAseq_ISE6etc_LibN228_N239.md # updated fastp (disregard "fastp.md")
```
### 2a. fastqc
```sh
fastqc.md
fastQC_SE_BGI_sRNAseq_ISE6etc_LibN228_N239.md # updated fastQC (disregard "fastqc.md")
```
### 3a multiQC
```sh
multiqc.md
multiqc_preprocess_SE_sRNAseq_N228toN239.md # updated multiQC (disregard "multiqc.md")
```

## II. Mapping of reads to Hlo and viral genomes

### 2a. Map reads to H. longicornis genome
```sh
STARmapping_SE_sRNAseqHlo_to_HaeL2018genome.md
```
### 2b. Map reads to viral genome
```sh
Bowtie_map_sRMAseqHlo_to_LGTVgenome.md
BowtieMap_SE_sRNAseqHlo_to_HAZVgenome_update.md
```
