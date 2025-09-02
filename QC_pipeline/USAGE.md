# MMC Analysis Pipeline Usage

## Overview

This guide describes how to use this directory to run the arcas pipeline, perform QC and filtering, and carry out Seurat processing to create finalised Seurat objects ready for downstream analysis.

Make sure you have first connected to the HPC before running anything - filepaths rely on this being the case! (ssh shs-sdf01)
---

## Input Files

The following information applies to directories started by copying the 'template/' base to conduct analysis

### FASTQs

Place your fastq.gz files in the 'data/reads/' folder

### 'samples.tsv' (in 'data/')

| Column     | Description                                                                 |
|------------|-----------------------------------------------------------------------------|
|'sample_id' | Unique sample identified (e.g. 'BLD001_cntl_R1'                             |
|'library_id | Typically 'MMC'                                                             |
|'well'      | Well sample was loaded in according to sample loading sheet (e.g. 'A1').    |
|            | Multiple wells should be comma separated without spaces (e.g. 'F3,F4,F5')   |

Example:

| sample_id     | library_id     | well     |
|---------------|----------------|----------|
|BLD001_cntl_R1 | MMC            | A1       |
|BLD001_LPS_R1  | MMC            | F2,F3,F4 |

### 'config/config.yaml'

- The user must edit the 'datadir' argument to the **full absolute path** of their analysis directory (the folder containing 'scripts/', 'data/', 'outputs/' and 'config/')
	- Example: 
		'''
		datadir: /mnt/odap-beegfs/a008/my_analysis_folder
		'''
- Other parameters are preconfigured and should not need modification

---

## Pipeline Scripts

The following information applies to directories started from the 'template/' base to conduct analysis

- **00.1_make_run_sheet.sh**
	Generates 'run_sheet.tsv' based on FASTQ files in 'data/reads/'. Must be run first.

- **00.2_dag.sh** (optional)
	Generates a DAG of the arcas pipeline (files saved in 'DAG/' which is created).
	Use after preparing all inputs (FASTQs, samples.tsv, run_sheet.tsv, config.yaml) to verify pipeline correctness.

- **01.1_run_pipeline.slurm**
	Runs the main arcas pipeline (long runtime: 2+ days).
	Requires 'samples.tsv', 'run_sheet.tsv', reads and config ready.

- ** 02.1_QC_and_filtering_pipeline.slurm**
	Performs quality control and filtering (runtime: several hours)
	Requires full outputs from 01.1
	
- **03.1_seurat_processing_pipeline.slurm**
	Performs Seurat downstream analysis (runtime: several hours)
	Requires initial Seurat list produced in 02.1

---

## SLURM Settings

- Users should update '--error' and '--output' file paths in all '.slurm' scripts to reflect their own directory (replace 'template' with the actual folder name).
- Adjust runtime and memory allocations as needed depending on dataset size.

---

## Usage Summary

1. Copy the 'template/' and rename it for your analysis

2. Add FASTQ files to 'data/reads/' and place your 'samples.tsv' in 'data/'

3. Edit 'config/config.yaml' to set the 'datadir' to your new analysis directory.

**For all next steps, edit the SLURM SETTINGS (--error --output and maybe --mem --time) and USER SETTINGS sections in the scripts before running**
**Assumed commands are being run from the analysis directory (that contains 'scripts/', 'data/' etc)

4. Run preprocessing script:
	'''
	sh scripts/00_preprocessing/00.1_make_run_sheet.sh
	'''
5. Optional: run the DAG script:
	'''
	sh scripts/00_preprocessing/00.2_dag.sh
	'''
6. Submit the main pipeline SLURM job:
	'''
	sbatch scripts/01_acras_pipeline/01.1_run_pipeline.slurm
	'''
	*** If this stops before the end of the run, to rerun you'll need to add --unlock to the very end of the command inside 01.1_run_pipeline.slurm, rerun, then remove this and rerun again
7. Once complete, submit QC and filtering:
	'''
	sbatch scripts/02_QC_and_filtering/02.1_QC_and_filtering_pipeline.slurm
	'''
8. Finally, submit Seurat processing:
	'''
	sbatch scripts/03_seurat_processing/03.1_seurat_processing_pipeline.slurm
	'''

---

## Multiple Runs Processing

If your analysis includes samples from **multiple sequencing runa (e.g. multiple flowcells)** and you want to process each individually (e.g. to compare flowcell efficiency), process each run individually through 0.1, then combine for downstream steps

Here's how:

1. Place only reads from the first run in 'data/reads/'
2. Add the corresponding 'samples.tsv' in 'data/'
3. Run '00.1', '00.2' and '01.1' as detailed above
4. Rename the files
	'''
	mv data/reads/ data/reads_runname1/
	mv data/samples.tsv data/samples_runname1.tsv
	mv data/run_sheet.tsv data/run_sheet_runname1.tsv
	'''
5. Add reads from the second run to a new 'data/reads/'
6. Add the second sample sheet as 'data/samples.tsv'
7. Update '01.1' USER SETTINGS with a new runname
8. Rerun '00.1' and '01.1'

This will result in:
	
outputs/arcas_pipeline/
|
|--runname1
|--runname2

Downstream steps will process these together, and the final metadata will include a 'flowcell' column indicating origin.

---

## Dependencies

- Pipelines use Singularity containers and conda environments pre-installed on the system
- If you encounter issues with dependencies, please contact the ODAP team for support.

---

*End of README*
