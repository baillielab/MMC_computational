rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")


library(Seurat)
library(SeuratDisk)

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

out_dirs<-paste(out_dirs, "/anndata_output")

seurat_list<-readRDS(paste(base, "/", out_dirs, "/intermediate/seurat_list_filtered_normalised_variable_doublet.rds", sep = ""))

# Create the output directory if it doesn't already exist
if (!dir.exists(paste(base, "/", outdir, sep = ""))) {
  dir.create(paste(base, "/", outdir, sep = ""))
}



### DEFINE FUNCTION ###

# Current SeuratObject version creates Assay5 - SeuratDisk cannot process this; needs back conversion
# Extract variable features and metadata info from current seurat object
# Extract count matrix and normalised matrix to create legacy assay
# Create .h5Seurat, convert to anndata (.h5ad)
convert_seurat_to_h5ad<-function(seu, filename_prefix) {
  var_features<-VariableFeatures(seu)
  meta_data<-seu@meta.data
  
  counts_mat<-GetAssayData(seu, assay = "RNA", slot = "counts")
  norm_data<-GetAssayData(seu, assay = "RNA", slot = "data")
  
  legacy_assay<-CreateAssayObject(counts = counts_mat)
  
  seu[["RNA"]]<-legacy_assay
  seu[["RNA"]]@data<-norm_data
  
  VariableFeatures(seu)<-var_features
  seu@meta.data<-meta_data
  
  h5seurat_file<-paste0(filename_prefix, ".h5Seurat")
  SaveH5Seurat(seu, filename = h5seurat_file, overwrite = TRUE)
  
  h5ad_file<-paste0(filename_prefix, ".h5ad")
  Convert(h5seurat_file, dest = "h5ad", overwrite = TRUE)
  
  message("Saved and converted: ", filename_prefix)
}



### PROCESS DATA ###

# For each object in seurat list, create anndata named with sample name
for (i in seq_along(seurat_list)) {
  seu<-seurat_list[[i]]
  sample_name<-names(seurat_list)[i]
  convert_seurat_to_h5ad(seu, filename_prefix = paste(base, "/", outdir, "/", sample_name, sep = ""))
}
