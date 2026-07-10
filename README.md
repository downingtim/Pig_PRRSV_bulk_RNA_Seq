# Pig_PRRSV_bulk_RNA_Seq

Analysis of longitudinal bulk RNA-Seq data from pigs experimentally infected with Porcine Reproductive and Respiratory Syndrome Virus (PRRSV). This repository contains the workflow used for RNA-Seq quality control, transcript quantification, differential expression analysis, and comparison of host transcriptional responses with PRRSV abundance.

---

## Repository Structure

```text
Pig_PRRSV_bulk_RNA_Seq/

├── FIGURES/
│   ├── PC1.PC2.pdf
│   └── PC3.PC4.pdf
│
├── TABLES/
│   └── DE_transcripts.csv
│
├── scripts/
│   ├── BSP462_PRRSV.Rmd
│   ├── Snakefile.kallisto
│   ├── Snakemake.qc
│   ├── abundanceplot.R
│   ├── indexscript.sh
│   ├── kraken.illumina.sh
│   └── kraken.illumina.extract.sh
│
└── README.md
```

---

## Overview

Pigs were sampled over the course of a PRRSV infection experiment and subjected to bulk RNA sequencing. The analysis workflow consisted of:

1. Read quality control and filtering
2. Transcriptome quantification with Kallisto
3. Transcript annotation with Ensembl BioMart
4. Differential expression analysis using Limma-Voom
5. Temporal modelling of transcript abundance
6. Identification of PRRSV reads using Kraken2
7. Correlation of pig transcript expression with viral abundance

---

## Experimental Design

Samples were collected at:

| Day post infection |
|-------------------|
| 0 |
| 3 |
| 7 |
| 10 |
| 14 |
| 21 |
| 28 |
| 35 |

Three pigs were sampled at each post-infection timepoint.

---

## Quality Control

Quality control and read trimming were performed using:

- FastQC
- MultiQC
- Fastp
- FASTX Toolkit

Workflow:

```bash
snakemake -s scripts/Snakemake.qc
```

The QC workflow:

- Generates FastQC reports for raw reads
- Performs quality filtering and adapter removal
- Generates FastQC reports for filtered reads
- Summarises all QC reports using MultiQC

---

## Transcriptome Index Generation

Transcriptome references were generated using `kb-python`.

Script:

```bash
scripts/indexscript.sh
```

Input:

- Ensembl FASTA
- Ensembl GTF

Included transcript classes:

- protein_coding
- lncRNA
- lincRNA
- antisense
- immunoglobulin genes
- T-cell receptor genes

Output:

```text
index.idx
t2g.txt
```

---

## Transcript Quantification

Transcript abundance estimation was performed using Kallisto.

Workflow:

```bash
snakemake -s scripts/Snakefile.kallisto
```

Kallisto outputs for each sample:

```text
abundance.tsv
abundance.h5
```

A simple abundance distribution plot can be generated using:

```bash
Rscript abundanceplot.R abundance.tsv abundance.png
```

---

## Principal Component Analysis

PCA was performed on transcript TPM values generated with Sleuth.

Outputs currently included in this repository:

```text
FIGURES/PC1.PC2.pdf
FIGURES/PC3.PC4.pdf
```

Timepoints were coloured using an 8-level Viridis palette corresponding to:

```text
D0
D3
D7
D10
D14
D21
D28
D35
```

---

## Differential Expression Analysis

Differential expression analysis was performed using:

- tximport
- edgeR
- limma
- voom

The design matrix treated each sampling day as a separate factor:

```r
design <- model.matrix(~0 + condition)
```

Pairwise comparisons against day 0 included:

```text
D3_vs_D0
D7_vs_D0
D10_vs_D0
D14_vs_D0
D21_vs_D0
D28_vs_D0
D35_vs_D0
```

Normalization and modelling workflow:

```r
filterByExpr()
calcNormFactors()
voom()
lmFit()
eBayes()
```

---

## Differentially Expressed Transcripts

Differentially expressed transcripts were defined as:

```text
FDR ≤ 0.05
|logFC| > 1.39
```

Output table:

```text
TABLES/DE_transcripts.csv
```

This table contains transcript-level differential expression statistics from the PRRSV time-course experiment.

---

## Transcript Annotation

Transcript annotations were retrieved from Ensembl using BioMart:

```r
biomaRt::getBM()
```

Annotations include:

- Ensembl transcript ID
- Transcript version
- Ensembl gene ID
- External gene symbol
- Transcript biotype
- Gene description

Reference archive:

```text
https://sep2025.archive.ensembl.org
```

---

## Temporal Modelling

Time-dependent expression changes were examined using natural splines:

```r
ns(day, df = 3)
```

Models were fitted using Limma:

```r
fitsp <- lmFit(v, model.matrix(~ ns(day, df=3)))
fitsp <- eBayes(fitsp, trend=TRUE)
```

Genes exhibiting significant temporal behaviour were identified using:

```r
topTableF()
```

---

## Viral Read Identification

PRRSV reads were identified using Kraken2.

Scripts:

```text
scripts/kraken.illumina.sh
scripts/kraken.illumina.extract.sh
```

PRRSV taxonomic ID:

```text
28344
```

Pig taxonomic ID:

```text
9822
```

KrakenTools was used to extract reads assigned to PRRSV:

```text
extract_kraken_reads.py
```

Outputs consist of FASTQ files containing reads classified as PRRSV.

---

## Host–Virus Correlation Analysis

Mean viral expression was calculated per sample and per timepoint from Kraken-derived PRRSV reads.

For each differentially expressed pig transcript:

1. Mean expression was calculated at each day.
2. Expression profiles were compared against the viral trajectory.
3. Pearson correlation coefficients were calculated.

This analysis was used to identify transcripts whose kinetics most closely followed PRRSV abundance over the infection time course.

---

## Major R Packages

```r
aggregation
biomaRt
dplyr
edgeR
ggplot2
ggpubr
ggrepel
ggvenn
ggVennDiagram
limma
readr
rhdf5
sleuth
tidyr
tximport
```

---

## Running the Analysis

### Quality Control

```bash
snakemake -s scripts/Snakemake.qc
```

### Transcript Quantification

```bash
snakemake -s scripts/Snakefile.kallisto
```

### Differential Expression and Downstream Analyses

```r
rmarkdown::render("scripts/BSP462_PRRSV.Rmd")
```

---

## Outputs

### Figures

```text
PC1.PC2.pdf
PC3.PC4.pdf
```

### Tables

```text
DE_transcripts.csv
```

---

## Contact

For questions regarding the workflow or analysis, please open a GitHub issue.

