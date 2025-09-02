# a008 - MMC Project Folder - README

## Overview

This folder contains all MMC-related analysis work. It includes various directories for raw data, reference materials, pipeline scripts, and outputs. 
Below is a description of the main folders and usage instructions.

The directory is intended for the processing if snRNA-seq data using the arcas-pipeline, followed by quality control (QC), filtering, and downstream Seurat-based analyses. 
The pipelines are intended to be modular, reproducible, and shared across users, with strict conventions to ensure consistency.

---

## Folder Structure and Descriptions

- 'arcas-main/'
	Contains the main arcas pipeline for processing FASTQ files into count matricies and downstream data.
	**Do not edit files here unless you know what you are doing**
	If edits are necessary, copy relevant files and edit the copy to avoid breaking the pipeline.

- 'template/'
	Base setup for running the arcas pipeline, QC and filtering, and Seurat steps.
	This folder should be copied and renamed when running new analyses.
	In their copied, renamed folder, users add reads to 'data/reads/', a 'samples.tsv' file to 'data/', and edit the 'config/config.yaml' file (only the 'datadir' argument) to point to the new analysis folder.
	Then run the pipeline scripts in order (follow instructions in the 'scripts' folder):
		- 'scripts/00_preprocessing/00.1_make_run_sheet.sh'
		- 'scripts/00_preprocessing/00.2_dag.sh' (optional but recommended)
		- 'scripts/01_arcas_pipeline/01.1_run_pipeline.slurm'
		- 'scripts/02_QC_and_filtering/02.1_QC_and_filtering_pipeline.slurm'
		- 'scripts/03_seurat_processing/03.1_seurat_processing_pipeline.slurm'

- 'pilot_pilot/'
	Analysis on the MMC pilot pilot data, completed following the steps outlined above.

- 'research/'
	Contains previous analyses, processed FASTQ files ('flowCellRuns/flowcell[12]/speciesSplit/'), and other research data.

- 'linked_data/'
	Metadata about sequencing runs, including sample loading table in '32890_Baillie_Kenneth/').

- 'scratch/' and user-specific directories
	For general directory setup and personal development; not required for general users.

- 'R_session_start.sh' 
	The script to start an R session on the HPC. Run sh R_session_start.sh then follow the instructions in the slurm file that appears.
---

## See USAGE.md for full details on how to run new analyses.

**End of README**
