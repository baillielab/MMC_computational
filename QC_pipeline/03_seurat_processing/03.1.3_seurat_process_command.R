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
library(future)
library(future.apply)

# Set up parallel plan
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

plan("multisession", workers = 16)
options(future.globals.maxSize = 250 * 1024^3)

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


seurat_list<-readRDS(paste(base, "/", out_dirs, "/intermediate/seurat_list_filtered_normalised_variable_doublet.rds", sep = ""))


### PROCESS DATA ###

# Separate tissues into separate objects
brain_samples<-grep("BRN", names(seurat_list), value = TRUE)
blood_samples<-grep("BLD", names(seurat_list), value = TRUE)
skin_samples<-grep("SKN", names(seurat_list), value = TRUE)

grouped_lists <- list(
  brain = seurat_list[brain_samples],
  blood = seurat_list[blood_samples],
  skin = seurat_list[skin_samples]
)

group_integrated <- list()

# To each list (tissue specific object) in list, process:
group_integrated <- future_lapply(names(grouped_lists), function(group_name) {
  group_list <- grouped_lists[[group_name]]
  
  # Downsample to reduce memory usage
  group_list<-lapply(group_list, function(x) {
    if (ncol(x) > 20000 ) {
      x<- subset(x, cells = sample(colnames(x), 20000))
    }
    return(x)
  })
  
  # Select features for integration
  features<-SelectIntegrationFeatures(object.list = group_list, nfeatures = 3000)
  
  # Scale and run PCA on each dataset to ensure shared features well represented before integration
  group_list <- lapply(group_list, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x<- RunPCA(x, features = features, verbose = FALSE)
    return(x)
  })
  
  # Find integration anchors (using reciprocal PCA for large datasets)
  anchors <- FindIntegrationAnchors(
    object.list = group_list,
    anchor.features = features,
    dims = 1:30,
    reduction = "rpca"
  )
  
  # Integrate data to single unified Seurat object
  integrated<- IntegrateData(anchorset = anchors, dims = 1:30, k.weight=80)
  integrated<- AddMetaData(integrated, metadata = group_name, col.name = "batch_group")
  
  return(integrated)
})

# Ensure names of elements in list are useful (matching input names)
names(group_integrated)<-names(grouped_lists)

# Save tissue specific integrated seurat objects
saveRDS(group_integrated$brain, paste(base, "/", out_dirs, "/intermediate/brain_integrated.rds", sep = ""))
saveRDS(group_integrated$blood, paste(base, "/", out_dirs, "/intermediate/blood_integrated.rds", sep = ""))
saveRDS(group_integrated$skin, paste(base, "/", out_dirs, "/intermediate/skin_integrated.rds", sep = ""))
