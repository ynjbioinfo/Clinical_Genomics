#!/bin/bash

# Set variables
REFERENCE_GENOME="hg38.fa"
REFERENCE_DICT="hg38.dict"
KNOWN_SITES_VCF="Homo_sapiens_assembly38.dbsnp138.vcf"
FASTQ_DIR="/path/to/fastq_bt"
OUTPUT_DIR="/path/to/output"
SAMPLE_NAME="demo"
FASTQ1="${FASTQ_DIR}/Demo_1.fq.gz"
FASTQ2="${FASTQ_DIR}/Demo_2.fq.gz"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Step 1: Data Preprocessing

# 1.1 FASTQC - Quality control of raw FASTQ files
echo "Running FASTQC..."
fastqc -o ${OUTPUT_DIR} ${FASTQ1} ${FASTQ2}

# 1.2 FASTP - Trimming adapter sequences and low-quality bases
echo "Running FASTP..."
fastp -i ${FASTQ1} -I ${FASTQ2} -o ${OUTPUT_DIR}/Demo_1_trimmed.fastq.gz -O ${OUTPUT_DIR}/Demo_2_trimmed.fastq.gz --detect_adapter_for_pe

# 1.3 Run FASTQC again after trimming
echo "Running FASTQC on trimmed files..."
fastqc -o ${OUTPUT_DIR} ${OUTPUT_DIR}/Demo_1_trimmed.fastq.gz ${OUTPUT_DIR}/Demo_2_trimmed.fastq.gz

# Step 2: Read Alignment

# 2.1 Index Reference Genome
echo "Indexing reference genome with BWA..."
bwa index ${REFERENCE_GENOME}

# 2.2 Create Reference Dictionary
echo "Creating reference dictionary with GATK..."
gatk CreateSequenceDictionary -R ${REFERENCE_GENOME} -O ${REFERENCE_DICT}

# 2.3 Align reads using BWA MEM
echo "Aligning reads with BWA MEM..."
bwa mem -t 4 -R "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tSM:${SAMPLE_NAME}" ${REFERENCE_GENOME} ${OUTPUT_DIR}/Demo_1_trimmed.fastq.gz ${OUTPUT_DIR}/Demo_2_trimmed.fastq.gz > ${OUTPUT_DIR}/${SAMPLE_NAME}.paired.sam

# Step 3: Post-Alignment Processing

# 3.1 Mark Duplicates and Sort
echo "Marking duplicates with GATK..."
gatk MarkDuplicatesSpark -I ${OUTPUT_DIR}/${SAMPLE_NAME}.paired.sam -O ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_dedup_reads.bam

# 3.2 Base Quality Score Recalibration (BQSR)
echo "Running BaseRecalibrator with GATK..."
gatk BaseRecalibrator -I ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_dedup_reads.bam -R ${REFERENCE_GENOME} --known-sites ${KNOWN_SITES_VCF} -O ${OUTPUT_DIR}/${SAMPLE_NAME}_recal_data.table

echo "Applying BQSR with GATK..."
gatk ApplyBQSR -I ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_dedup_reads.bam -R ${REFERENCE_GENOME} --bqsr-recal-file ${OUTPUT_DIR}/${SAMPLE_NAME}_recal_data.table -O ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam

# Step 4: Variant Calling

# 4.1 Call Variants
echo "Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller -R ${REFERENCE_GENOME} -I ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_dedup_bqsr_reads.bam -O ${OUTPUT_DIR}/${SAMPLE_NAME}_raw_variants.vcf

# 4.2 Filter Variants
echo "Filtering SNPs and INDELs with GATK..."
gatk SelectVariants -R ${REFERENCE_GENOME} -V ${OUTPUT_DIR}/${SAMPLE_NAME}_raw_variants.vcf --select-type SNP -O ${OUTPUT_DIR}/${SAMPLE_NAME}_raw_snps.vcf
gatk SelectVariants -R ${REFERENCE_GENOME} -V ${OUTPUT_DIR}/${SAMPLE_NAME}_raw_variants.vcf --select-type INDEL -O ${OUTPUT_DIR}/${SAMPLE_NAME}_raw_indels.vcf

# Apply filters
echo "Applying variant filters with GATK..."
gatk VariantFiltration \
    -R ${REFERENCE_GENOME} \
    -V ${OUTPUT_DIR}/${SAMPLE_NAME}_raw_snps.vcf \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_snps.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

# Select variants that pass filters
echo "Selecting variants that pass filters with GATK..."
gatk SelectVariants \
    --exclude-filtered \
    -V ${OUTPUT_DIR}/${SAMPLE_NAME}_filtered_snps.vcf \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_analysis_ready_snps.vcf

# Step 5: Annotation

# 5.1 Convert VCF to ANNOVAR format
echo "Converting VCF to ANNOVAR format..."
docker run -it --rm -v ${OUTPUT_DIR}:/data bioinfochrustrasbourg/annovar:latest perl ./convert2annovar.pl -format vcf4 /data/${SAMPLE_NAME}_filtered_snps.vcf -outfile /data/${SAMPLE_NAME}_filtered_snps.avinput

# 5.2 Download required annotation files
echo "Downloading ANNOVAR reference files..."
wget -P ${OUTPUT_DIR} http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz
wget -P ${OUTPUT_DIR} http://www.openbioinformatics.org/annovar/download/hg38_refGeneMrna.fa.gz
wget -P ${OUTPUT_DIR} http://www.openbioinformatics.org/annovar/download/hg38_refGeneVersion.txt.gz

# Move reference files into ANNOVAR directory
echo "Preparing ANNOVAR database..."
mkdir -p ${OUTPUT_DIR}/annovar/humandb
mv ${OUTPUT_DIR}/hg38_refGene.txt.gz ${OUTPUT_DIR}/annovar/humandb/
mv ${OUTPUT_DIR}/hg38_refGeneMrna.fa.gz ${OUTPUT_DIR}/annovar/humandb/
mv ${OUTPUT_DIR}/hg38_refGeneVersion.txt.gz ${OUTPUT_DIR}/annovar/humandb/

echo "Pipeline completed successfully!"
