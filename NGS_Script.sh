

echo "------Welcome in NGS pipeline We wil start step by step-------"
echo "Prepreation files- Download the file e.g. Here i am taking two pair end  RCL2341_1.fq.gz  and RCL2341_2.fq.gz"

echo "download the fastq files"
echo "Reference files -wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
echo "gunzip hg38.fa.gz"

echo "# index ref - .fai file before running haplotype caller"
echo "samtools faidx hg38.fa #Done"


echo "# ref dict - .dict file before running haplotype caller"
echo "gatk CreateSequenceDictionary R=hg38.fa O=hg38.dict #Done"
#for running the docker 
sudo docker run broadinstitute/gatk:latest

# for checking the commands 
sudo docker run broadinstitute/gatk:latest gatk HaplotypeCaller --help


# download known sites files for BQSR from GATK resource bundle
echo "wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf #Done"
echo "wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx #Done"

# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc RCL2341_1.fq.gz  #Done
fastqc RCL2341_2.fq.gz  #Done

# No trimming required, quality looks okay.


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

echo "# BWA index reference "
echo "bwa index hg38.fa #Done"


# BWA alignment
bwa mem -t 4 -R "@RG\tID:RCL2341\tPL:ILLUMINA\tSM:RCL2341" hg38.fa RCL2341_1.fq.gz RCL2341_2.fq.gz > RCL2341.paired.sam #Done



# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I RCL2341.paired.sam -O RCL2341_sorted_dedup_reads.bam  #Done

# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I RCL2341_sorted_dedup_reads.bam -R hg38.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf  -O recal_data.table #done


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I RCL2341_sorted_dedup_reads.bam -R hg38.fa  --bqsr-recal-file recal_data.table -O RCL2341_sorted_dedup_bqsr_reads.bam #Done

# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=hg38.fa I=RCL2341_sorted_dedup_bqsr_reads.bam O=alignment_metrics.txt #Done
gatk CollectInsertSizeMetrics INPUT=RCL2341_sorted_dedup_bqsr_reads.bam OUTPUT=insert_size_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf #Not Done R is required for this


# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R hg38.fa -I RCL2341_sorted_dedup_bqsr_reads.bam -O raw_variants.vcf #Done



# extract SNPs & INDELS

gatk SelectVariants -R hg38.fa -V raw_variants.vcf --select-type SNP -O raw_snps.vcf #Done
gatk SelectVariants -R hg38.fa -V raw_variants.vcf --select-type INDEL -O raw_indels.vcf #Done


# Filter SNPs    #Done
gatk VariantFiltration \
	-R hg38.fa  \
	-V raw_snps.vcf \
	-O filtered_snps.vcf \
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

# ----------------------------------------------
# STEP 6: Call Variants - Variant annotation
# ----------------------------------------------
