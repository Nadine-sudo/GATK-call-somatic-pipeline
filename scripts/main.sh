#!/bin/bash

# Main control script for GATK somatic variant calling pipeline

# Load common functions and variables
source $(dirname "$0")/common.sh

# Parse command-line arguments
while getopts "1:2:d:t:r:h" opt
do
	case $opt in
		"1") CONTROL1=$OPTARG;;
		"2") CONTROL2=$OPTARG;;
		"d") WORK_DIR=$OPTARG;;
		"t") TUMOR=$OPTARG;;
	    "r") REF_SPECIES=$OPTARG;;
		"h") usage; exit 0;;
		"?") echo "Invalid option: -$OPTARG" >&2; usage; exit 1;;
		":") echo "Option -$OPTARG requires an argument." >&2; usage; exit 1;;  
	esac
done

WORK_DIR=${WORK_DIR:-$(pwd)}
SAMPLE_NAME=$(basename $TUMOR)
DATA=$(dirname "$CONTROL1")
CONTROL1_NAME=$(basename $CONTROL1)
CONTROL2_NAME=$(basename $CONTROL2)
REF_SPECIES=${REF_SPECIES:-"human"}

# Check if required parameters are provided
if [[ -z "$CONTROL1" || -z "$TUMOR" ]]; then
	echo "Error: Tumor and at least one control sample are required." >&2
	usage
	exit 1
fi

# Select reference genome based on species
if [[ "$REF_SPECIES" == "rat" ]]; then
  export REFERENCE=$REFERENCE_RAT
  export RESOURCE=$RESOURCE_RAT
else
  export REFERENCE=$REFERENCE_HUMAN
  export RESOURCE=$RESOURCE_HUMAN
fi
export REF_SPECIES
export WORK_DIR

# Add error handling for missing scripts
if [[ ! -f "$(dirname "$0")/preprocess.sh" || ! -f "$(dirname "$0")/variant_calling.sh" ]]; then
    echo "Error: Required scripts (preprocess.sh or variant_calling.sh) are missing." >&2
    exit 1
fi

# Main pipeline logic
if [ ! -e "${WORK_DIR}/${CONTROL1_NAME}/6.BQSR/${CONTROL1_NAME}.sorted.markdup.BQSR.bam" ]; then
    if [[ -n "$CONTROL2_NAME" ]]; then
        SAMPLE_LIST="${SAMPLE_NAME} ${CONTROL1_NAME} ${CONTROL2_NAME}"
    else
        SAMPLE_LIST="${SAMPLE_NAME} ${CONTROL1_NAME}"
    fi
    for SAMPLE in $SAMPLE_LIST
    do
        mkdir -p ${WORK_DIR}/${SAMPLE}
        cd ${WORK_DIR}/${SAMPLE}
        $(dirname "$0")/preprocess.sh ${DATA} ${SAMPLE}
    done
    if [[ -n "$CONTROL2" ]]; then
        $(dirname "$0")/variant_calling.sh ${SAMPLE_NAME} ${CONTROL1_NAME} ${CONTROL2_NAME}
    else
        $(dirname "$0")/variant_calling.sh ${SAMPLE_NAME} ${CONTROL1_NAME} ""
    fi
else
    cd ${WORK_DIR}/${SAMPLE_NAME}
    if [[ -n "$CONTROL2" ]]; then
        $(dirname "$0")/variant_calling.sh ${SAMPLE_NAME} ${CONTROL1_NAME} ${CONTROL2_NAME}
    else
        $(dirname "$0")/variant_calling.sh ${SAMPLE_NAME} ${CONTROL1_NAME} ""
    fi
fi
