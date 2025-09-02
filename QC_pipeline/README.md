# Scripts Directory

This folder contains all scripts required to run the MMC snRNA-seq processing pipeline.

## Pipeline Overview

Run scripts in this order:

1. '00.1_make_run_sheet.sh' - Generates 'run_sheet.tsv'
2. '00.2_dag.sh' - (Optional) Runs a DAG to check pipeline will execute correctly
3. '01.1_run_pipeline.slurm' - Runs the full arcas pipeline to produce count matricies etc
4. '02.1_QC_and_filtering.pipeline.slurm' - Performs QC including mapping statistics, busfile QC metrics, and initial Seurat QC and filtering
5. '03.1_seurat_processing_pipeline.slurm' - Performs downstream Seurat-based processing to generate final tissue specific Seurat objects ready for analysis

---

## What to Edit

For each script:

- Update the **USER SETTINGS** section
	Insstructions on what to update and how are in the script itself
- Update the **SLURM SETTINGS** section
	Update --error and --output log file paths (change 'template' to your actual directory name).
	Update --time and --mem as needed

---

## Do Not Edit

- Files like '01.1.1_Snakefile_bustools', '02.1.1_...' etc are downstream scripts called by the main pipelines
	These should **not** be edited unless you are debugging or modifying the workflow itself.
