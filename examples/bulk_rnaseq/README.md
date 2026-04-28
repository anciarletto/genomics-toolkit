# Bulk RNA-Seq Pipeline

This folder contains a reproducible pipeline for bulk RNA-seq analysis using publicly available human data (GEO). It covers all steps from raw FASTQ files through alignment, quantification, differential expression, visualization, and pathway enrichment analysis.

---

## Features

- End-to-end RNA-seq pipeline
- Parallelized HISAT2 alignment
- Automated reference genome and annotation download
- Handles paired-end data
- Reproducible and modular design
- Tested with GEO dataset: GSE326480 (breast cancer brain metastasis)

---

## Pipeline Overview

The workflow includes:

1. Download raw sequencing data and metadata from SRA.
2. Quality control with FastQC and MultiQC.
3. Optional adapter trimming using Cutadapt.
4. Alignment to the human genome (GRCh38) using HISAT2.
5. Quantification of gene expression using featureCounts.
6. Differential expression analysis with DESeq2.
7. Visualization of results (volcano plots, heatmaps, PCA, etc.).
8. Conduct pathway enrichment analysis with GSEA

FASTQ  
&nbsp;&nbsp;  ↓   
FastQC → MultiQC  
&nbsp;&nbsp;  ↓  
(Optional) Cutadapt  
&nbsp;&nbsp;  ↓  
HISAT2 alignment  
&nbsp;&nbsp;  ↓  
Sorted BAM files  
&nbsp;&nbsp;  ↓  
featureCounts  
&nbsp;&nbsp;  ↓  
Gene count matrix  
&nbsp;&nbsp;  ↓  
DESeq2  
&nbsp;&nbsp;  ↓  
Differential expression + plots  
&nbsp;&nbsp;  ↓  
GSEA   


---

## Scripts and Descriptions

- `download_fastq.sh` – Downloads sequencing data from SRA and converts to FASTQ.
- `get_metadata.R` - Generates a DESeq2-compatible sample metadata table.  
- `run_fastqc.sh` – Performs quality control on FASTQ files.  
- `run_cutadapt.sh` *(optional)* – Trims adapters and low-quality bases.  
- `run_hisat2.sh` – Aligns reads to the GRCh38 reference genome and outputs sorted BAM files.  
- `run_featurecounts.sh` – Quantifies gene-level counts from BAM files.  
- `run_deseq2.R` – Performs differential expression analysis using DESeq2, pathway enrichment analysis with GSEA, and generates visualizations.

---

## Usage

Run the pipeline step by step:

### 1. Download FASTQ files and Metadata
./scripts/download_fastq.sh  
./scripts/get_metadata.R

### 2. Quality control
./scripts/run_fastqc.sh

### 3. Summarize QC
multiqc qc_reports/ -o multiqc_reports/

### 4. (Optional) Trim reads
./scripts/run_cutadapt.sh

### 5. Align reads to reference genome
./scripts/run_hisat2.sh

### 6. Generate gene counts
./scripts/run_featurecounts.sh

### 7. Differential expression analysis
Rscript scripts/run_deseq2.R

---

## Example dataset

This pipeline was tested using GEO dataset:
GSE326480 (breast cancer brain metastasis RNA-seq)

---

## Requirements

### Core Tools
- HISAT2 – read alignment
- samtools – BAM processing
- Subread package (featureCounts) – gene-level quantification
- FastQC – read quality control
- MultiQC – aggregate QC reports

### Data Retrieval
- SRA Toolkit (`prefetch`, `fasterq-dump`) – download sequencing data
- pigz (or gzip) – FASTQ compression

### R Environment
- R (≥ 4.0 recommended)
- Bioconductor packages:
  - DESeq2
  - fgsea
- CRAN packages:
  - msigdbr
  - ggplot2
  - dplyr
  - readr
  - pheatmap
  - RColorBrewer
  - stringr

### Other
- bash shell (for running pipeline scripts)

### Reference Files
- Reference genome (FASTA)
- Gene annotation (GTF)

### Input Data
- SRR accession list (used for downloading FASTQ files)
- Sample metadata file (experimental conditions, replicates)

### Gene Sets
- MSigDB Hallmark gene sets (accessed via msigdbr R package)

---

## Notes

- Data is paired-end
- Library is unstranded (`-s 0` in featureCounts)
- Reference genome: GRCh38
