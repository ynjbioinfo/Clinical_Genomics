# Clinical_Genomics
Clinical_Genomics_WES_Pipeline
# Genomic Data Analysis Pipeline

## Introduction

This repository provides a comprehensive pipeline for genomic data analysis, specifically focusing on whole-genome sequencing and variant calling using the Genome Analysis Toolkit (GATK). It includes steps for data preprocessing, read alignment, post-alignment processing, variant calling, and annotation. This pipeline utilizes Docker for GATK to ensure consistent and reproducible results. 

The provided commands and scripts guide users through the entire workflow from raw FASTQ files to annotated variant calls, using widely adopted bioinformatics tools and best practices.

## Requirements

To run this pipeline, ensure the following are installed:

- **Docker**: For running the GATK container and other tools in isolated environments.
- **GATK**: The Genome Analysis Toolkit for variant discovery in high-throughput sequencing data.
- **BWA**: For aligning short DNA sequences to a reference genome.
- 
- **SAMtools**: For manipulating SAM/BAM files.
- **FASTQC**: For quality control of raw sequencing data.
- **FASTP**: For quality trimming and adapter removal from FASTQ files.
- **ANNOVAR**: For functional annotation of genetic variants.

## Tools and Docker Setup

### Docker Setup

**Running the GATK Docker Container**

To run the GATK Docker container, use the following command:

```bash
sudo docker run broadinstitute/gatk:latest

To view the available commands in the GATK container:


sudo docker run broadinstitute/gatk:latest gatk HaplotypeCaller --help

Author
This pipeline and documentation were developed by Yogesh Joshi. For questions or feedback, please contact Yogesh Joshi.

Pipeline Overview
1. Data Preprocessing
Objective: Prepare raw sequencing data for alignment by performing quality control and trimming.

FASTQC:

Tool: FASTQC
Result Interpretation: Provides a quality report of raw FASTQ files, identifying potential issues such as low quality scores or adapter contamination.
FASTP:

Tool: FASTP
Result Interpretation: Trims adapter sequences and low-quality bases from the reads. Outputs trimmed FASTQ files ready for alignment.
2. Read Alignment
Objective: Align the preprocessed reads to a reference genome.

BWA:
Tool: BWA (Burrows-Wheeler Aligner)
Result Interpretation: Aligns short DNA sequences to the reference genome. Outputs a SAM file, which contains the aligned read information.
3. Post-Alignment Processing
Objective: Process the alignment files to prepare them for variant calling.

MarkDuplicatesSpark:

Tool: GATK MarkDuplicatesSpark
Result Interpretation: Identifies and marks duplicate reads in the SAM/BAM file, which helps to reduce biases in variant calling. Outputs a sorted, deduplicated BAM file.
Base Quality Score Recalibration (BQSR):

Tool: GATK BaseRecalibrator and ApplyBQSR
Result Interpretation: Recalibrates base quality scores based on known variant sites to correct systematic errors. Outputs a recalibrated BAM file with improved accuracy for variant calling.
4. Variant Calling
Objective: Identify genetic variants from the aligned reads.

HaplotypeCaller:

Tool: GATK HaplotypeCaller
Result Interpretation: Calls variants (SNPs and INDELs) from the aligned BAM files. Outputs a VCF (Variant Call Format) file with raw variant calls.
SelectVariants:

Tool: GATK SelectVariants
Result Interpretation: Filters and separates variants into SNPs and INDELs. Outputs filtered VCF files for further analysis.
VariantFiltration:

Tool: GATK VariantFiltration
Result Interpretation: Applies filters to remove low-confidence variants based on specific quality metrics. Outputs a filtered VCF file with high-confidence variants.
5. Annotation
Objective: Annotate the identified variants to understand their biological significance.

ANNOVAR:
Tool: ANNOVAR
Result Interpretation: Annotates variants with information about gene function, disease associations, and other biological features. Outputs annotated variant files with functional and clinical annotations.

