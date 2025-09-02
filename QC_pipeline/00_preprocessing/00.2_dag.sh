#!/bin/bash 

# === USER SETTINGS ===

# Set base paths

# File path to where analysis will take place, same as datadir in config file
BASE=/mnt/odap-beegfs/a008/pilot_pilot_new

# Snakefile options for end orientation. 
# For paired end (automatic orientation detection) write: Snakefile_paired
# For unpaired reads, or uncertainty with orientation write: Snakefile_unpaired
SNAKEFILE_TYPE=Snakefile_unpaired

# =================================================



# Set base paths - invariable
ROOT=/mnt/odap-beegfs/a008
CONT=/mnt/odap-beegfs/software/singularity-images

# Set specific paths
CONT_ARCAS=${CONT}/arcas-container-v2.sif
CONT_GVIZ=${CONT}/graphviz-container.sif
SNAKEFILE=${ROOT}/arcas-main/workflow/${SNAKEFILE_TYPE}
OUTPUT=${BASE}/DAG
CONFIG=${BASE}/config/config.yaml

# === CREATE DAG ===

# Create DAG directory (if doesn't exist)
mkdir -p "${OUTPUT}"

# Create dot file
singularity exec --bind ${ROOT},${CONT} ${CONT_ARCAS} snakemake --dag -s ${SNAKEFILE} --configfile ${CONFIG} --cores 16 --use-conda --use-singularity -d ${OUTPUT} > ${OUTPUT}/dag.dot

# dot to pdf
singularity exec --bind ${ROOT},${CONT} ${CONT_GVIZ} dot -Tpdf ${OUTPUT}/dag.dot -o ${OUTPUT}/dag.pdf

# dot to png
singularity exec --bind ${ROOT},${CONT} ${CONT_GVIZ} dot -Tpng ${OUTPUT}/dag.dot -o ${OUTPUT}/dag.png

