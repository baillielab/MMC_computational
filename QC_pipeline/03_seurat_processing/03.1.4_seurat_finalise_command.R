rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(Seurat)
library(ggplot2)
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

# Read in tissue split integrated seurat objects from previous step
brain<-readRDS(paste(base, "/", out_dirs, "/intermediate/brain_integrated.rds", sep = ""))
blood<-readRDS(paste(base, "/", out_dirs, "/intermediate/blood_integrated.rds", sep = ""))
skin<-readRDS(paste(base, "/", out_dirs, "/intermediate/skin_integrated.rds", sep = ""))


### DEFINE FUNCTIONS ###
# Set up functions as each process needs to be done for all tissue types

# Add tissue type, condition (LPS or cntl) and sample.id to metadata of seurat object
# Return the updated object
process_metadata <- function(seurat_obj, tissue_type) {
  seurat_obj$tissue<-tissue_type

  seurat_obj$condition <-ifelse(
    grepl("cntl", seurat_obj$orig.ident, ignore.case = TRUE),
    "control",
    "LPS"
  )

  seurat_obj$sample_id<-seurat_obj$orig.ident

  return(seurat_obj)

}

# Set assay to integrated, scale data, run PCA, run UMAP, find neighbours and find clusters
# Return updated seurat object
dim_red<-function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "integrated"

  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

  return(seurat_obj)
}

# Create UMAP plots grouped by seurat clusters, condition, flowcell, and sample
# Create PCA plot of condition
# Create violin plots of QC metrics (nFeature_RNA, nCount_RNA, percent.mt)
# Save all of these
plots<-function(seurat_obj, outstem) {
  clusters<-DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")
  dir.create(paste(base, "/", out_dirs, "/figures/", outstem, sep = ""))
  ggsave(paste(base, "/", out_dirs, "/figures/", outstem, "/", outstem, "_umap_clusters.png", sep = ""), clusters, bg = "white", width = 10, height=6)

  condition<-DimPlot(seurat_obj, reduction = "umap", group.by = "condition")
  ggsave(paste(base, "/", out_dirs, "/figures/", outstem, "/", outstem, "_umap_condition.png", sep = ""), condition, bg = "white", width = 10, height=6)

  flowcell<-DimPlot(seurat_obj, reduction = "umap", group.by = "flowcell")
  ggsave(paste(base, "/", out_dirs, "/figures/", outstem, "/", outstem, "_umap_flowcell.png", sep = ""), flowcell, bg = "white", width = 10, height=6)

  sample<-DimPlot(seurat_obj, reduction = "umap", group.by = "sample_id")
  ggsave(paste(base, "/", out_dirs, "/figures/", outstem, "/", outstem, "_umap_sample.png", sep = ""), sample, bg = "white", width = 10, height=6)

  PCA<-DimPlot(seurat_obj, reduction = "pca", group.by = "condition")
  ggsave(paste(base, "/", out_dirs, "/figures/", outstem, "/", outstem, "_PCA_clusters.png", sep = ""), PCA, bg = "white", width = 10, height=6)

  vln<-VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", pt.size = 0.01, ncol = 3)
  ggsave(paste(base, "/", out_dirs, "/figures/", outstem, "/", outstem, "_VlnPlot.png", sep = ""), vln, bg = "white", width = 10, height=6)

}

# Run Idents on object
# Find all markers (save rds of markers)
# Extract metadata (save csv)
# Return updated seurat object
clusters<-function(seurat_obj, outstem) {
  seurat_obj$cluster<-Idents(seurat_obj)
  markers<-FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(markers, paste(base, "/", out_dirs, "/metadata/", outstem, "_all_markers.rds", sep = ""))
  metadata<-seurat_obj@meta.data
  write.csv(metadata, paste(base, "/", out_dirs, "/metadata/", outstem, "_seurat_metadata.csv", sep = ""), row.names = TRUE)

  return(seurat_obj)
}



### PROCESS THE DATA ###

# Run each of the functions on each of the tissue specific seurat objects

brain<-process_metadata(brain, "brain")
blood<-process_metadata(blood, "blood")
skin<-process_metadata(skin, "skin")

brain<-dim_red(brain)
blood<-dim_red(blood)
skin<-dim_red(skin)

plots(brain, "brain")
plots(blood, "blood")
plots(skin, "skin")

brain<-clusters(brain, "brain")
blood<-clusters(blood, "blood")
skin<-clusters(skin, "skin")


# Save the updated tissue specific seurat objects
# These have now been fully QC'ed and processed, ready for analysis
saveRDS(brain, paste(base, "/", out_dirs, "/analysis_ready/brain_integrated_final.rds", sep = ""))
saveRDS(blood, paste(base, "/", out_dirs, "/analysis_ready/blood_integrated_final.rds", sep = ""))
saveRDS(skin, paste(base, "/", out_dirs, "/analysis_ready/skin_integrated_final.rds", sep = ""))


