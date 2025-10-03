#!/bin/bash

# Common functions and variables for GATK somatic variant calling pipeline

function usage(){
    NAME=$(basename $0)
    cat << EOF
Usage:
${NAME}[OPTIONS]      

OPTIONS:
        
-1  input sample path of control 1
-2  input sample path of control 2
-d  working directory(default:  .)
-t  path of tumor sample
-r  reference species (human or rat, default: human)
-h  display this help message

EOF
}

# Load environments and set paths 
if ! source path/to/activate path/to/gatk; then
    echo "Error: Failed to activate GATK environment." >&2
    exit 1
fi

BWA="path/to/bwa"
GATK="python3 path/to/gatk"
MANTA="path/to/manta"
STRELKA="path/to/strelka"
DATABASE="path/to/database"
TEMP="path/to/temp"
VEP="path/to/vep"
FASTP="path/to/fastp"
SAMTOOLS="path/to/samtools"

# Add more reference genome paths as needed
REFERENCE_RAT="path/to/reference/rat/rn7/rn7.fa"
REFERENCE_HUMAN="path/to/reference/human/hg38/hg38.fa"

# Public resources for GATK BQSR
RESOURCE_HUMAN="path/to/public/data/GATK/known-sites/hg38"
RESOURCE_RAT="path/to/public/data/GATK/known-sites/rat"
