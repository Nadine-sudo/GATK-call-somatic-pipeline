#!/bin/bash

# Preprocess raw data to generate BAM files

# Load common functions and variables
source $(dirname "$0")/common.sh

# Add error handling for missing tools
if ! command -v fastp &> /dev/null || ! command -v samtools &> /dev/null; then
    echo "Error: Required tools (fastp or samtools) are not installed or not in PATH." >&2
    exit 1
fi

# Ensure output directories exist
mkdir -p 1.fastp 2.bwa 3.samtool 4.flagstat 5.markdup 6.BQSR

function preprocess(){
    # $1: data directory
    # $2: sample name
    local DATA_DIR=$1
    local SAMPLE=$2

    # Only look for raw data
    local RAW1=$(ls ${DATA_DIR}/${SAMPLE}/*_1.fq.gz 2>/dev/null)
    local RAW2=$(ls ${DATA_DIR}/${SAMPLE}/*_2.fq.gz 2>/dev/null)

    if [[ -n "$RAW1" && -n "$RAW2" ]]; then
    mkdir -p 1.fastp
    fastp -w 8 \
        -i $RAW1 \
        -I $RAW2 \
        -o 1.fastp/$(basename ${RAW1}) \
        -O 1.fastp/$(basename ${RAW2}) \
        --html 1.fastp/${SAMPLE}.html \
        --json 1.fastp/${SAMPLE}.json \
        1> 1.fastp/fastp.log 2>&1

    # Mapping to reference genome
    mkdir -p 2.bwa
    $BWA mem -t 12 -M -Y \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        ${REFERENCE} \
        1.fastp/$(basename ${RAW1}) \
        1.fastp/$(basename ${RAW2}) \
        -o 2.bwa/${SAMPLE}.sam \
        1> 2.bwa/bwa.log 2>&1

    # Convert SAM to BAM
    mkdir -p 3.samtool
    samtools view \
        -bS 2.bwa/${SAMPLE}.sam \
        -o 3.samtool/${SAMPLE}_2.bam \
        1> 3.samtool/samtool.log 2>&1

    # Merge BAMs if 3.samtool/${SAMPLE}_1.bam exists, otherwise remove
    if [ -f "3.samtool/${SAMPLE}_1.bam" ]; then
        samtools merge -@ 8 -h 3.samtool/${SAMPLE}_2.bam \
        3.samtool/${SAMPLE}.bam \
        3.samtool/${SAMPLE}_1.bam \
        3.samtool/${SAMPLE}_2.bam
    else
        mv 3.samtool/${SAMPLE}_2.bam 3.samtool/${SAMPLE}.bam
    fi
    samtools sort -@ 8 3.samtool/${SAMPLE}.bam -o 3.samtool/${SAMPLE}.sorted.bam 1> 3.samtool/sort.log 2>&1

    # BAM QC statistics
    mkdir -p 4.flagstat
    samtools flagstat 3.samtool/${SAMPLE}.sorted.bam > 4.flagstat/${SAMPLE}.sorted.stat

    # Mark duplicates and create BAM index
    mkdir -p 5.markdup
    $GATK MarkDuplicates \
        -I 3.samtool/${SAMPLE}.sorted.bam \
        -O 5.markdup/${SAMPLE}.sorted.markdup.bam \
        -M 5.markdup/${SAMPLE}.sorted.markdup_metrics.txt \
        --CREATE_INDEX true \
        --TMP_DIR $TEMP \
        1> 5.markdup/markdup.log 2>&1

        # Base Quality Score Recalibration (BQSR)
        mkdir -p 6.BQSR

        # Set known-sites arguments based on reference species
        if [[ "$REF_SPECIES" == "rat" ]]; then
            KNOWN_SITES=(
                "--known-sites ${RESOURCE}/HRDP_smp146_HPJoint_GATK4_rnBN7_INDELs_HF_PASS_snpEff.vcf.gz"
                "--known-sites ${RESOURCE}/HRDP_smp146_HPJoint_GATK4_rnBN7_SNPs_HF_PASS_snpEff.vcf.gz"
            )
        else
            KNOWN_SITES=(
                "--known-sites ${RESOURCE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
                "--known-sites ${RESOURCE}/Homo_sapiens_assembly38.known_indels.vcf.gz"
                "--known-sites ${RESOURCE}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
                "--known-sites ${RESOURCE}/Homo_sapiens_assembly38.dbsnp138.vcf"   
            )
        fi
        
        $GATK BaseRecalibrator \
            -R ${REFERENCE} \
            -I 5.markdup/${SAMPLE}.sorted.markdup.bam \
            ${KNOWN_SITES[@]} \
            -O 6.BQSR/recal_data_${SAMPLE}.table \
            1> 6.BQSR/BaseRecalibrator.log 2>&1

        $GATK ApplyBQSR \
            -R ${REFERENCE} \
            --bqsr-recal-file 6.BQSR/recal_data_${SAMPLE}.table \
            -I 5.markdup/${SAMPLE}.sorted.markdup.bam \
            -O 6.BQSR/${SAMPLE}.sorted.markdup.BQSR.bam \
            1> 6.BQSR/ApplyBQSR.log 2>&1
    else
        echo "No valid raw FASTQ files found for sample $SAMPLE in $DATA_DIR" >&2
        return 1
    fi
}

# Call the preprocess function with arguments
preprocess $1 $2
