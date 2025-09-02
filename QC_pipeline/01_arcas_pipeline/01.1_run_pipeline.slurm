#!/bin/bash 

# === SLURM SETTINGS ===

#SBATCH --job-name=arcas-pipeline
#SBATCH --time=96:00:00
#SBATCH --mem=100G
#SBATCH --output=/mnt/odap-beegfs/a008/pilot_pilot_new/logs/%x/%x_%j.out
#SBATCH --error=/mnt/odap-beegfs/a008/pilot_pilot_new/logs/%x/%x_%j.err
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32

# =======================



# === USER SETTINGS ===

# Set base paths

# File path to where analysis will take place, same as datadir in config file 
BASE=/mnt/odap-beegfs/a008/pilot_pilot_new
# The run name - this will be the name of your output folder
RUNNAME=FC1_unpaired

# Snakefile options for end orientation.                                        
# For paired end (automatic orientation detection) write: Snakefile_paired             
# For unpaired reads, or uncertainty with orientation write: Snakefile_unpaired
SNAKEFILE_TYPE=Snakefile_unpaired

# ==========================



# Set base paths - invariable
ROOT=/mnt/odap-beegfs/a008
CONT=/mnt/odap-beegfs/software/singularity-images

#Set specific paths
CONT_ARCAS=${CONT}/arcas-container-v2.sif
SNAKEFILE=${ROOT}/arcas-main/workflow/${SNAKEFILE_TYPE}
OUTPUT=${BASE}/outputs/arcas_pipeline/${RUNNAME}
CONFIG=${BASE}/config/config.yaml


# === RUN PIPELINE ===
mkdir -p ${BASE}/outputs/arcas_pipeline
singularity exec --bind ${ROOT},${CONT} ${CONT_ARCAS} snakemake -s ${SNAKEFILE} --configfile ${CONFIG} --cores 64 --use-conda --use-singularity --rerun-incomplete -d ${OUTPUT}
