rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(scDblFinder)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

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

seurat_list<-readRDS(paste(base, "/", out_dirs, "/intermediate/seurat_list_filtered_normalised_variable.rds", sep = ""))

### DEFINE FUNCTION ###

# Convert seurat object to SingleCellExperiment, run scDblFinder and add outputs 
# class and score to seurat object. Output the object
run_scdblfinder<-function(seurat_obj) {
  
  sce<-as.SingleCellExperiment(seurat_obj)
  
  set.seed(123)
  sce<-scDblFinder(sce)
  
  seurat_obj$scDblFinder_class<-sce$scDblFinder.class
  seurat_obj$scDblFinder_score<-sce$scDblFinder.score

  return(seurat_obj)
  
}


### PROCESS DATA ###

# Run the function for each object in the seurat list
seurat_list<-lapply(seurat_list, run_scdblfinder)


### ANALYSE DATA ###

# Check the percentage of doublets for each sample
sapply(seurat_list, function(seu) {
  round(prop.table(table(seu$scDblFinder_class))["doublet"] * 100, 2)
})

# Create data frame identifying singlets/doublets for each element of seurat_list
df<-do.call(rbind, lapply(names(seurat_list), function(nm) {
  data.frame(
    sample = nm,
    status = seurat_list[[nm]]$scDblFinder_class
  )
}))

# Create a barblot showing the percentage of each sample that are singlets and doublets
# Save it
bars<-ggplot(df, aes(x = sample, fill = status)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion of cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)); bars

ggsave(paste(base, "/", out_dirs, "/figures/qc/doublets.png", sep = ""),
       bars, bg="white", width = 10, height = 6)

# Create a violin plot of confidence (score) for each cell as singlet / doublet
# Save it
violin<-VlnPlot(seurat_list[[1]], features = "scDblFinder_score", group.by = "scDblFinder_class")

ggsave(paste(base, "/", out_dirs, "/figures/qc/doublet_confidence.png", sep = ""),
       violin, bg="white", width = 10, height = 6)


# Remove doublets from each object
seurat_list<-lapply(seurat_list, function(seu) {
  subset(seu, subset = scDblFinder_class == "singlet")
})

# Save the updated seurat_list
saveRDS(seurat_list, 
        paste(base, "/", out_dirs, "/intermediate/seurat_list_filtered_normalised_variable_doublet.rds", sep = ""))

