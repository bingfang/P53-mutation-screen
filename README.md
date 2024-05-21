---
title: "P53 mutation screen"
author: "Bingfang Xu"
date: "04/22/2024"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Run snakemake pipeline

```
# dry run
snakemake --dry-run check_target
snakemake -np check_target

# run
snakemake check_target --cores 1
# snakemake -np

# unlock
# snakemake --unlock

# analyze workflow
#snakemake --lint

```

# Install trimmomatic

```
#test conda 
conda env list
conda update conda

#add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#install trimmamatic
conda install -c bioconda trimmomatic
```

# Install biopython

```
conda install biopython

```

# Trimmatic commend lines

```
#paired end trimming WITHOUT trimming start seqence
#input files: *.fastq.gz or *.fastq
java -jar /$MyPATH/trimmomatic-0.39-2/trimmomatic.jar PE -phred33 P53-samples_S1_L001_R1_001.fastq P53-samples_S1_L001_R2_001.fastq P53-samples_S1_L001_R1_trimmed.fastq P53-samples_S1_L001_R1_trimmed_unpair.fastq P53-samples_S1_L001_R2_trimmed.fastq P53-samples_S1_L001_R2_trimmed_unpair.fastq TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

```

# Subset the fastq files (optional)

```
# if reads per samples is more than 100000, subset the fastq file to reduce computation time.
/$MyPATH2/seqtk sample -s100 P53_samples_R1_trimmed.fastq 500000 > P53_samples_R1_sub_trimmed.fastq
/$MyPATH2/seqtk sample -s100 P53_samples_R2_trimmed.fastq 500000 > P53_samples_R2_sub_trimmed.fastq

```

# Run python scripts one by one

Debarcode, genotyping, and mutation rate calculation

```
# cd to the script dir
cd $MyPATH3/src

# input sample barcode info, then run debarcode script.
python3 _0debarcode_Miseq.py

# check the filter for mutation rate, then run genotyping script.
python3 _1_P53_genotyped_mouseTrp53.py

# Open Trp53 snapGene, and align sequences from "Check_seq_all.txt" 
# input ~10bp WT, mutated sequences, and run the script of mutation rate calculation.
python3 _2CheckSNP_mouseTrp53.py
```