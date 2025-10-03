# GATK Somatic Mutation Analysis Pipeline  
A modular, automated workflow for somatic SNV detection for tumor-normal paired samples or tumor-only samples, fully aligned with **GATK Best Practices for Somatic Mutation Calling**. 
This pipeline streamlines end-to-end analysis—from raw FastQ data to functionally annotated variants—with reusable scripts for scalability and reproducibility.


## 1. Pipeline Core Modules  
The workflow is split into 4 key modules, each managed by a dedicated script (stored in `scripts/`). Below is the module dependency and data flow:  

```mermaid
flowchart TD
    A[Raw FastQ<br>(Tumor + Normal)] -->|Input| B[preprocess.sh<br>Data Preprocessing]
    B -->|Output: Recalibrated BAMs| C[variant_calling.sh<br>Somatic Calling]
    C -->|Output: Filtered VCF| D[intersect_and_annotate.sh<br>Variant Post-Analysis]
    E[common.sh<br>Shared Utilities] --> B
    E --> C
    E --> D
    F[main.sh<br>Pipeline Entry] --> B
    F --> C
    F --> D
```

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
