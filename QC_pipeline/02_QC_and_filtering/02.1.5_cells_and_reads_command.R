rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(data.table)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)
library(jsonlite)

# Set file paths
base=args[1]
in_dirs=args[2]
subdir<-args[3]
files<-Sys.glob(paste(base, "/", in_dirs, sep = ""))
in_dirs<-NA
for (i in 1:length(files)) {
  in_dirs<-c(in_dirs, tail(strsplit(files[i], "/")[[1]], n=1))
}
in_dirs<-in_dirs[-1]
in_dirs<-paste("outputs/arcas_pipeline/", in_dirs, "/", subdir, sep = "")
out_dirs=args[4]

# Sample names - import from run and sample sheet
sample_names<-read.table(args[5])

# Read in list of seurat objects from previous step
seurat_list<-readRDS(paste(base, "/outputs/seurat_processing/intermediate/seurat_list_filtered_normalised_variable.rds", sep = ""))

# Fix sample names
sample_names<-sample_names[,1]
sample_names<-sample_names[-1]
if (length(grep("Saccular", sample_names)) != 0) {
  sample_names<-sample_names[-c(grep("Saccular", sample_names))] # Removing saccular as not part of MMC analysis
}



### PROCESS DATA ###

# Create empty dataframe to populate
reads_per_cell<-data.frame(sample = NA,
                           flowcell = NA,
                           reads = NA,
                           cells = NA,
                           rpc = NA)


for (j in 1:length(in_dirs)) {
  # Loop over all input directories
  fc<-in_dirs[j]
  
  print(paste("Working on flowcell", j))
  
  for (i in 1:length(sample_names)) {
    # Loop over all samples
    sample<-sample_names[i]
    print(paste("Processing sample", sample, "number", i, "of", length(sample_names)))
    
    # Add flowcell to name to distinguish samples with same name
    name<-paste(sample, "_FC", j, sep = "")
    
    # Read sample specific json
    json<-read_json(paste(base, "/", fc, "/", sample, "/kallisto.run_info.json", sep = ""))
    
    # Extract info from json
    total_reads<-json$n_processed
    
    # Extract number of cells in sample (post filtering) from seurat list
    cells<-ncol(seurat_list[[name]])
    
    # Calculate number of reads per cell
    calc<-round(total_reads/cells)
    
    # Fill data frame
    add<-data.frame(sample = sample,
                    flowcell = j,
                    reads = total_reads,
                    cells = cells,
                    rpc = calc)
    
    # Add to existing data frame
    reads_per_cell<-rbind(reads_per_cell, add)
    
  }
  
}

# Remove empty row from setting up dataframe
reads_per_cell<-reads_per_cell[-1,]

# Add tissue type
reads_per_cell$type<-NA

reads_per_cell$type[grep("BRN", reads_per_cell$sample)]<-"Brain"
reads_per_cell$type[grep("BLD", reads_per_cell$sample)]<-"Blood"
reads_per_cell$type[grep("SKN", reads_per_cell$sample)]<-"Skin"

# Overwrite scientific notation option
options(scipen=999)

# Plot number of reads per cell per sample, split by flowcell
# Save it
rpc<-ggplot(reads_per_cell, aes(x=sample, y=rpc, fill = factor(flowcell))) +
  geom_col(position = "dodge")+
  labs(title = "Reads per cell", y = "Reads per Cell", x = "Sample", fill = "Flowcell") +
  theme_minimal() + 
  coord_cartesian(ylim=c(0,50000))+
  geom_hline(yintercept=20000)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); rpc
ggsave(paste(base, "/", out_dirs, "/reads_per_cell.png", sep = ""), rpc, bg = "white", width = 10, height=6)

# Boxplot to compare number of reads per cell between flowcells
reads_fc<-ggplot(reads_per_cell, aes(x = factor(flowcell), y = rpc)) +
  geom_boxplot()+
  geom_jitter(width=0.2, alpha = 0.6) +
  coord_cartesian(ylim=c(0,50000))+
  labs(title = "Reads per cell per Flowcell", y = "Cells", x = "Flowcell") +
  theme_minimal(); reads_fc
ggsave(paste(base, "/", out_dirs, "/reads_per_cell_fc.png", sep = ""), reads_fc, bg = "white", width = 10, height=6)

# Boxplot to compare number of reads per cell between tissues
reads_tissue<-ggplot(reads_per_cell, aes(x = type, y = rpc)) +
  geom_boxplot()+
  geom_jitter(width=0.2, alpha = 0.6) +
  coord_cartesian(ylim=c(0,50000))+
  labs(title = "Reads per cell per tieeus", y = "Cells", x = "Tissue") +
  theme_minimal(); reads_tissue
ggsave(paste(base, "/", out_dirs, "/reads_per_cell_tissue.png", sep = ""), reads_tissue, bg = "white", width = 10, height=6)


# Plot of number of cells per sample, split by flowcell
cells<-ggplot(reads_per_cell, aes(x=sample, y=cells, fill = factor(flowcell))) +
  geom_col(position = "dodge")+
  labs(title = "Cells per sample", y = "Cells", x = "Sample", fill = "Flowcell") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); cells
ggsave(paste(base, "/", out_dirs, "/cells_per_sample.png", sep = ""), cells, bg = "white", width = 10, height=6)

# Boxplot to compare number of cells per sample across flowcells
cells_fc<-ggplot(reads_per_cell, aes(x = factor(flowcell), y = cells)) +
  geom_boxplot()+
  geom_jitter(width=0.2, alpha = 0.6) +
  labs(title = "Cells per sample per Flowcell", y = "Cells", x = "Flowcell") +
  theme_minimal(); cells_fc
ggsave(paste(base, "/", out_dirs, "/cells_per_sample_fc.png", sep = ""), cells_fc, bg = "white", width = 10, height=6)

# Boxplot to compare number of cells per sample across tissues
cells_tissue<-ggplot(reads_per_cell, aes(x = type, y = cells)) +
  geom_boxplot()+
  geom_jitter(width=0.2, alpha = 0.6) +
  labs(title = "Cells per sample per tissue", y = "Cells", x = "Tissue") +
  theme_minimal() ; cells_tissue
ggsave(paste(base, "/", out_dirs, "/cells_per_sample_tissue.png", sep = ""), cells_tissue, bg = "white", width = 10, height=6)


