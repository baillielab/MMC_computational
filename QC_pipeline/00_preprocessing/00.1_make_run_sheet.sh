#!/bin/bash

# === USER SETTINGS ===

# Set base file path (where you will run the pipeline from - the same as datadir for the config file)
BASE=/mnt/odap-beegfs/a008/pilot_pilot_new

# Fastq pattern
## Define how to extract the sample ID from the filename:
## Character that splits different parts of the filename (e.g. for test_sample_test_R1.fastq.gz this is _ so SPLIT_CHAR="_")
SPLIT_CHAR="_"

## Identify parts to remove / keep from the start and end of the filename EXCLUDING R[12].fastq.gz
## (e.g. for test_sample_test_R1.fastq.gz we want to remove 1 (test) from the start, and 1 (test) from the end so REMOVE_START=1 and REMOVE_END=1)
REMOVE_START=1
REMOVE_END=2

# ============================================================



# Set specific file paths
READS=${BASE}/data/reads    
OUTPUT=${BASE}/data/run_sheet.tsv

# Write header
echo -e "run_id\tlibrary_id\tfastq1\tfastq2" > ${OUTPUT}

# Declare output maos
declare -A fastq1_map
declare -A fastq2_map

library_id=MMC

# Loop through all R1/R2 files in READS
for file in ${READS}/*_R[12].fastq.gz; do
	filename=$(basename "$file")
	
	# === CLEAN AND PARSE FILENAME ===

	# Strip _R[12].fastq.gz ending
	name_no_ext=$(echo "$filename" | sed -E 's/_R[12]\.fastq\.gz$//')
	# Collapse any double underscores
	name_no_ext="${name_no_ext//__/_}"
	
	# Split using specified character
	split_safe="${name_no_ext//${SPLIT_CHAR}/ }"
	read -ra parts <<< "$split_safe"

	# Compute part boundaries
	total_parts=${#parts[@]}
	start=$REMOVE_START
	end=$((total_parts - $REMOVE_END))
	
	# Extract run_id parts
	run_id_parts=("${parts[@]:start:end-start}")
	run_id=$(IFS="${SPLIT_CHAR}"; echo "${run_id_parts[*]}")
	
	# Map R1 and R2 files to run_id
	read_type=$(echo "$filename" | grep -o 'R[12]')

	if [[ "$read_type" == "R1" ]]; then
		fastq1_map["$run_id"]="$file"
	elif [[ "$read_type" == "R2" ]]; then
		fastq2_map["$run_id"]="$file"
	fi
done

# === OUTPUT TSV ===

# fastq file paths should be from the directory the pipeline will be run in - not the full file path - remove
RELATIVE=$BASE/

for run_id in "${!fastq1_map[@]}"; do
	fastq1="${fastq1_map[$run_id]}"
	fastq2="${fastq2_map[$run_id]}"

	# Strip leading path (make relative)
	fastq1_rel="${fastq1#$RELATIVE}"
	fastq2_rel="${fastq2#$RELATIVE}"

	echo -e "${run_id}\t${library_id}\t${fastq1_rel}\t${fastq2_rel}" >> ${OUTPUT}
done

echo "Output written to: $OUTPUT"

