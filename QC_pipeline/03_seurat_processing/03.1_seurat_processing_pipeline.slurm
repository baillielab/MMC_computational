#!/bin/bash                                                                     

# === SLURM SETTINGS ===
#SBATCH --job-name=seurat_processing                                                      
#SBATCH --time=72:00:00                                                         
#SBATCH --mem=300G                                                               
#SBATCH --output=/mnt/odap-beegfs/a008/pilot_pilot_new/logs/%x/%x_%j.out           
#SBATCH --error=/mnt/odap-beegfs/a008/pilot_pilot_new/logs/%x/%x_%j.err            
#SBATCH --ntasks=1                                                              
#SBATCH --cpus-per-task=32  

# =========================


# === USER SETTINGS ===

# Set base paths                                                                

# File path to where analysis will take place, same as datadir in config file
BASE=/mnt/odap-beegfs/a008/pilot_pilot_new



# =================================================

# Set base paths                                                                
CONT=/mnt/odap-beegfs/software/singularity-images

#Set specific paths                                                             
CONT_R=${CONT}/mmc-rstudio-container-v3.sif

# Set config arguments
SUBDIR=MMC
IN_DIR=outputs/arcas_pipeline/* 
OUTDIR=outputs/seurat_processing

# === PROCESS SEURAT OBJECTS  ===

# Identify and remove doublets
singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/03_seurat_processing/03.1.1_doublets.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}"

# Create anndata
singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/03_seurat_processing/03.1.2_anndata_creation.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}"

# Integrate into tissue specific objects
singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/03_seurat_processing/03.1.3_seurat_process_command.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}"

# Extract metadata, process integrated objects (scale, PCA, UMAP, clusters), plot QC, PCA and UMAP, find all markers
# Save analysis ready seurat objects for each tissue
mkdir -p "${BASE}/${OUTDIR}/metadata"
mkdir -p "${BASE}/${OUTDIR}/analysis_ready"
singularity exec --bind ${BASE},${CONT} ${CONT_R} Rscript ${BASE}/scripts/03_seurat_processing/03.1.4_seurat_finalise_command.R "${BASE}" "${IN_DIR}" "${SUBDIR}" "${OUTDIR}"
