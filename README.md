# GATK Somatic Mutation Analysis Pipeline  
A modular, automated workflow for somatic SNV detection for tumor-normal paired samples or tumor-only samples, fully aligned with **GATK Best Practices for Somatic Mutation Calling**. 
This pipeline streamlines end-to-end analysis—from raw FastQ data to functionally annotated variants—with reusable scripts for scalability and reproducibility.


## 1. Pipeline Core Modules  
The workflow is split into 4 key modules, each managed by a dedicated script (stored in `scripts/`). 
The pipeline follows a sequential workflow with clear module dependencies:
- Raw FastQ files (tumor and normal samples) are first processed by preprocess.sh to generate recalibrated BAM files.
- These BAM files are then used as input for variant_calling.sh to perform somatic variant calling and produce filtered VCF files.
- The filtered VCFs are passed to intersect_and_annotate.sh for variant intersection, functional annotation, and final result generation.
All three core modules (preprocess.sh, variant_calling.sh, intersect_and_annotate.sh) rely on common.sh for shared utility functions (e.g., logging, file validation). The entire workflow is orchestrated by main.sh, which serves as the pipeline entry point and coordinates execution of all modules.

## 2. Repository File Structure
All files are organized for clarity, ensuring easy navigation and usage:
```
gatk-somatic-mutation-pipeline/
├── scripts/                     # Core executable scripts (run via main.sh or individually)
│   ├── common.sh                # Shared functions: logging, file validation, directory creation
│   ├── preprocess.sh            # Step 1: FastQC → BWA alignment → MarkDuplicates → BaseRecalibration
│   ├── variant_calling.sh       # Step 2: Mutect2 → FilterMutectCalls → Tumor purity estimation
│   ├── intersect_and_annotate.sh# Step 3: Variant intersection → VEP/ANNOVAR annotation → Stats
│   └── main.sh                  # Entry point: orchestrates all modules with user parameters
├── docs/                        # Documentation for setup, usage, and troubleshooting
│   └── environment.yaml           # Software dependencies (Conda/non-Conda) + reference genome prep
└──  README.md                    # Repository homepage (core instructions)
```

## 3. Quick Start Guide
### 3.1 Prerequisites
- Conda (Miniconda3/Anaconda3, recommended for dependency management)
- Git (for repository cloning)
### 3.2 Install Dependencies
Use the pre-configured Conda environment to avoid version conflicts:
```bash
# 1. Clone the repository
git clone https://github.com/your-username/gatk-somatic-mutation-pipeline.git
cd gatk-somatic-mutation-pipeline

# 2. Create and activate Conda environment (from docs/environment.yaml)
conda env create -f docs/environment.yaml
conda activate gatk-somatic
```
### 3.3 Data Preparation
To run the pipeline, you'll need the following data files. Here are recommended sources:
- Reference Genome (FASTA + Indexes)：
 - Before running the pipeline, you need to download the appropriate reference genome for your species. Reference genome must be indexed(using samtools faidx, bwa index and gatk CreateSequenceDictionary) before.
- Known-sites for BQSR：
 - Human Reference (hg38 example), use [GATK Resource Bundle - Homo_sapiens_assembly38.fasta](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38)
 - Other species (rat example), download from other database [RGD](https://download.rgd.mcw.edu/data_release/)
- Annotation Resources：
  - For human samples: Use GATK-compatible resources for Funcotator
  - For non-human species: Use VEP cache files matching your reference genome
### 3.4 Run the Pipeline
Execute the full workflow via main.sh with your input parameters:
```
bash scripts/main.sh \
  --tumor_fq1 /path/to/tumor_R1.fastq.gz \
  --tumor_fq2 /path/to/tumor_R2.fastq.gz \
  --normal_fq1 /path/to/normal_R1.fastq.gz \
  --normal_fq2 /path/to/normal_R2.fastq.gz \
  --ref /path/to/hg38.fa \
  --out_dir ./somatic_results \
  --threads 12
  ```
