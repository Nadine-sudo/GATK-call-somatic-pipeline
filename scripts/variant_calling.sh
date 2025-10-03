#!/bin/bash

# Variant calling script for GATK somatic variant calling pipeline

# Load common functions and variables
source $(dirname "$0")/common.sh

# Add error handling for missing tools
if ! command -v $GATK &> /dev/null; then
    echo "Error: GATK is not installed or not in PATH." >&2
    exit 1
fi

function mutect(){
    local SAMPLE_NAME=$1
    local CONTROL=$2

    mkdir -p 7.mutect/${SAMPLE_NAME}_${CONTROL}

    # Run Mutect2
    $GATK Mutect2 \
        -R ${REFERENCE} \
        -I 6.BQSR/${SAMPLE_NAME}.sorted.markdup.BQSR.bam \
        -I ../${CONTROL}/6.BQSR/${CONTROL}.sorted.markdup.BQSR.bam \
        -normal ${CONTROL} \
        --f1r2-tar-gz 7.mutect/${SAMPLE_NAME}_${CONTROL}/f1r2.tar.gz \
        -O 7.mutect/${SAMPLE_NAME}_${CONTROL}/${SAMPLE_NAME}_${CONTROL}.unfiltered.vcf.gz \
        --tmp-dir $TEMP \
        1> 7.mutect/${SAMPLE_NAME}_${CONTROL}/mutect.log 2>&1

    # Learn orientation bias artifacts
    $GATK LearnReadOrientationModel \
        -I 7.mutect/${SAMPLE_NAME}_${CONTROL}/f1r2.tar.gz \
        -O 7.mutect/${SAMPLE_NAME}_${CONTROL}/read-orientation-model.tar.gz \
        --tmp-dir $TEMP \
        1> 7.mutect/${SAMPLE_NAME}_${CONTROL}/orientation.log 2>&1

    # Set common germline variant sites arguments based on reference species
        if [[ "$REF_SPECIES" == "rat" ]]; then
            SITES=(
                "-V ${RESOURCE}/HRDP_48smp_HPJoint_gatk4_rnBN7_SNPs_HF_PASS.vcf.gz"
	            "-L ${RESOURCE}/rat.interval_list"
            )
        else
            SITES=(
                "-V ${RESOURCE}/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
	            "-L ${RESOURCE}/wgs_calling_regions.hg38.interval_list"
            )
        fi

	# summarize known germline variant resource
	for file in ${1} ${2}
	do
	    if [ -e "${WORK_DIR}/${file}/7.mutect/${file}.pileups.table" ]; then
	        continue
	    else
	        mkdir -p ${WORK_DIR}/${file}/7.mutect
	        $GATK GetPileupSummaries \
	        -I ${WORK_DIR}/${file}/6.BQSR/${file}.sorted.markdup.BQSR.bam \
	        ${SITES[@]} \
	        -O ${WORK_DIR}/${file}/7.mutect/${file}.pileups.table \
	        --tmp-dir $TEMP \
            1> ${WORK_DIR}/${file}/7.mutect/pileups.log 2>&1
	    fi
	done

	# estimate cross-sample contamination
	$GATK CalculateContamination \
	-I 7.mutect/${1}.pileups.table \
	-matched ../${2}/7.mutect/${2}.pileups.table \
	-O 7.mutect/${1}_${2}/contamination.table \
	--tmp-dir $TEMP \
    1> 7.mutect/${SAMPLE_NAME}_${CONTROL}/contamination.log 2>&1

    # Filter Mutect calls
    $GATK FilterMutectCalls \
        -R $REFERENCE \
        -V 7.mutect/${SAMPLE_NAME}_${CONTROL}/${SAMPLE_NAME}_${CONTROL}.unfiltered.vcf.gz \
        --contamination-table 7.mutect/${SAMPLE_NAME}_${CONTROL}/contamination.table \
        -ob-priors 7.mutect/${SAMPLE_NAME}_${CONTROL}/read-orientation-model.tar.gz \
        -O 7.mutect/${SAMPLE_NAME}_${CONTROL}/${SAMPLE_NAME}_${CONTROL}.filtered.vcf.gz \
        --tmp-dir $TEMP 

}

function strelka(){
    local SAMPLE_NAME=$1
    local CONTROL=$2

    mkdir -p 8.manta/${SAMPLE_NAME}_${CONTROL}

    # Run Manta for structural variant calling
    python2 ${MANTA}/configManta.py \
        --normalBam ../${CONTROL}/6.BQSR/${CONTROL}.sorted.markdup.BQSR.bam \
        --tumorBam 6.BQSR/${SAMPLE_NAME}.sorted.markdup.BQSR.bam \
        --referenceFasta $REFERENCE \
        --runDir 8.manta/${SAMPLE_NAME}_${CONTROL} \
        1> 8.manta/${SAMPLE_NAME}_${CONTROL}/config.log 2>&1

    python2 8.manta/${SAMPLE_NAME}_${CONTROL}/runWorkflow.py -j 8 \
        1> 8.manta/${SAMPLE_NAME}_${CONTROL}/execute.log 2>&1

    # Run Strelka for somatic variant calling
    mkdir -p 9.strelka/${SAMPLE_NAME}_${CONTROL}
    python2 ${STRELKA}/configureStrelkaSomaticWorkflow.py \
        --normalBam ../${CONTROL}/6.BQSR/${CONTROL}.sorted.markdup.BQSR.bam \
        --tumorBam 6.BQSR/${SAMPLE_NAME}.sorted.markdup.BQSR.bam \
        --referenceFasta ${REFERENCE} \
        --runDir 9.strelka/${SAMPLE_NAME}_${CONTROL} \
        --indelCandidates 8.manta/${SAMPLE_NAME}_${CONTROL}/results/variants/candidateSmallIndels.vcf.gz \
        1> 9.strelka/${SAMPLE_NAME}_${CONTROL}/config.log 2>&1

    python2 9.strelka/${SAMPLE_NAME}_${CONTROL}/runWorkflow.py -m local -j 16 \
        1> 9.strelka/${SAMPLE_NAME}_${CONTROL}/execute.log 2>&1
}

function call(){
    cd ${WORK_DIR}/${1}

    # Run Mutect and Strelka for CONTROL1
    mutect ${1} ${2}
    strelka ${1} ${2}

    # Check if CONTROL2 is provided
    if [[ -n "${3}" ]]; then
        mutect ${1} ${3}
        strelka ${1} ${3}
    fi

    
    # Call intersect and annotate script
    "$(dirname "$0")/intersect_and_annotate.sh" "$1" "$2" "$3"
    
}

# Call the call function with arguments
call $1 $2 $3
