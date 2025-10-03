#!/bin/bash

# Intersect and annotate variants

# Load common functions and variables
source $(dirname "$0")/common.sh

# Add error handling for missing tools
if ! command -v $GATK &> /dev/null || ! command -v $VEP/vcf2maf.pl &> /dev/null; then
    echo "Error: Required tools (GATK or vep) are not installed or not in PATH." >&2
    exit 1
fi


function intersect(){
    local SAMPLE_NAME=$1
    local CONTROL1=$2
    local CONTROL2=$3

    # Ensure required files exist for CONTROL1
    if [[ ! -f "7.mutect/${SAMPLE_NAME}_${CONTROL1}/${SAMPLE_NAME}_${CONTROL1}.filtered.vcf.gz" ]]; then
        echo "Error: Required Mutect VCF file for CONTROL1 is missing." >&2
        exit 1
    fi

    if [[ ! -f "9.strelka/${SAMPLE_NAME}_${CONTROL1}/results/variants/somatic.snvs.vcf.gz" ]]; then
        echo "Error: Required Strelka VCF file for CONTROL1 is missing." >&2
        exit 1
    fi

    # Intersect Mutect variants
    zgrep '^#' 7.mutect/${SAMPLE_NAME}_${CONTROL1}/${SAMPLE_NAME}_${CONTROL1}.filtered.vcf.gz > 7.mutect/${SAMPLE_NAME}.mutect.vcf


    if [[ -n "$CONTROL2" ]]; then
        if [[ ! -f "7.mutect/${SAMPLE_NAME}_${CONTROL2}/${SAMPLE_NAME}_${CONTROL2}.filtered.vcf.gz" ]]; then
            echo "Error: Required Mutect VCF file for CONTROL2 is missing." >&2
            exit 1
        fi

        awk 'NR==FNR{a[$1"\t"$2];next}(($1"\t"$2)in a)' \
            <(zcat 7.mutect/${SAMPLE_NAME}_${CONTROL1}/${SAMPLE_NAME}_${CONTROL1}.filtered.vcf.gz | grep 'PASS') \
            <(zcat 7.mutect/${SAMPLE_NAME}_${CONTROL2}/${SAMPLE_NAME}_${CONTROL2}.filtered.vcf.gz | grep 'PASS') >> 7.mutect/${SAMPLE_NAME}.mutect.vcf
    else
        zcat 7.mutect/${SAMPLE_NAME}_${CONTROL1}/${SAMPLE_NAME}_${CONTROL1}.filtered.vcf.gz | grep 'PASS' >> 7.mutect/${SAMPLE_NAME}.mutect.vcf
    fi

    # Intersect Strelka variants

    if [[ -n "$CONTROL2" ]]; then
        if [[ ! -f "9.strelka/${SAMPLE_NAME}_${CONTROL2}/results/variants/somatic.snvs.vcf.gz" ]]; then
            echo "Error: Required Strelka VCF file for CONTROL2 is missing." >&2
            exit 1
        fi

        zcat 9.strelka/${SAMPLE_NAME}_${CONTROL1}/results/variants/*gz | grep 'PASS' | sort -V > 9.strelka/${SAMPLE_NAME}_${CONTROL1}/results/variants/${SAMPLE_NAME}_${CONTROL1}_filter.vcf
        zcat 9.strelka/${SAMPLE_NAME}_${CONTROL2}/results/variants/*gz | grep 'PASS' | sort -V > 9.strelka/${SAMPLE_NAME}_${CONTROL2}/results/variants/${SAMPLE_NAME}_${CONTROL2}_filter.vcf
        
        awk 'NR==FNR {A[$1"\t"$2]; next}(($1"\t"$2) in A)' \
            9.strelka/${SAMPLE_NAME}_${CONTROL1}/results/variants/${SAMPLE_NAME}_${CONTROL1}_filter.vcf \
            9.strelka/${SAMPLE_NAME}_${CONTROL2}/results/variants/${SAMPLE_NAME}_${CONTROL2}_filter.vcf > 9.strelka/${SAMPLE_NAME}.strelka.vcf
    else
        zcat 9.strelka/${SAMPLE_NAME}_${CONTROL1}/results/variants/*gz | grep 'PASS' > 9.strelka/${SAMPLE_NAME}.strelka.vcf
    fi

    # Combine Mutect and Strelka intersections
    grep '^#' 7.mutect/${SAMPLE_NAME}.mutect.vcf > ${SAMPLE_NAME}.vcf
    awk 'NR==FNR{a[$1"\t"$2];next}(($1"\t"$2) in a)' \
        9.strelka/${SAMPLE_NAME}.strelka.vcf 7.mutect/${SAMPLE_NAME}.mutect.vcf >> ${SAMPLE_NAME}.vcf
}

function annotate(){
    if [[ "$REF_SPECIES" == "human" ]]; then
        $GATK Funcotator \
        -R $REFERENCE \
        --output-file-format MAF \
        --ref-version hg38 \
        --data-sources-path $DATABASE/funcotator_dataSources.v1.7.20200521s \
        -V ${1}.vcf \
        -O ${1}.maf
    else
        conda deactivate
        conda activate vep
        perl $VEP/vcf2maf.pl --input-vcf ${1}.vcf --output-maf ${1}.maf \
        --vep-path $VEP \
        --ref-fasta $REFERENCE \
        --vep-data $DATABASE/vep_database/rat/rn7 \
        --species rattus_norvegicus \
        --ncbi-build mRatBN7.2 \
        --cache-version 105 
    fi
}

# Call intersect and annotate functions
intersect $1 $2 $3
annotate $1
