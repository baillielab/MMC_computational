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
gene_names<-read.csv(paste(base, "/data/reference/t2g.csv", sep = ""))

# Fix sample names
sample_names<-sample_names[,1]
sample_names<-sample_names[-1]
if (length(grep("Saccular", sample_names)) != 0) {
  sample_names<-sample_names[-c(grep("Saccular", sample_names))] # Removing saccular as not part of MMC analysis
}


# ALL count.genes.names.txt FILES IDENTICAL - ONLY DO ONCE
# Read in gene ids to add to matrix
genes<-read.delim(paste(base, "/", in_dirs[1], "/", sample_names[1], "/count.genes.names.txt", sep = ""), header = FALSE)
genes$name<-genes$V1

# To use gene names instead of IDs:

# for (i in 1:length(genes$V1)) {
#   if (length(which(gene_names$gene == genes$V1[i])) != 0) {
#     genes$name[i]<-unique(gene_names$gene_name[which(gene_names$gene == genes$V1[i])])
#   }
# }


### FUNCTIONS ###
# Takes the count matrix count.genes.mature.mtx, and turns it into an annotated count matrix using gene ids from count.genes.names.txt and gene names from gene_names input table as
# row (gene) names and extracts barcodes from busfile to use as column (barcode) names
make_annotated_matrix<-function(flowcell, sample, gene_names) {
  # Read in count matrix (mature), annotate gene names
  counts<-readMM(paste(base, "/", flowcell, "/", sample, "/count.genes.mature.mtx", sep = ""))
  counts<-t(counts)
  
  rownames(counts)<-make.unique(genes$name)

  # Read in busfile to extract barcodes
  bus<-fread(paste(base, "/", flowcell, "/", sample, "/sorted_bus.txt", sep = ""))
  colnames(bus)<-c("barcode", "UMI", "ec", "count")
  
  barcodes<-unique(bus$barcode)
  
  if(length(barcodes) != ncol(counts)) {
    print("Number of barcodes does not match number of columns in count matrix!")
  } else {
    colnames(counts)<-barcodes
  }
  
  barcode_file<-file.path(paste(base, "/", flowcell, "/", sample, "/QC/barcode_filter.csv", sep = ""))
  
  if(file.exists(barcode_file)) {
    message("Using filtered barcodes for ", sample)
    
    filtered_barcodes<-fread(barcode_file, header = FALSE)$V1
    counts<-counts[,colnames(counts) %in% filtered_barcodes]
  }
  return(counts)
}

# Takes the count matrix and creates a seurat object, calculates % mitochondrial genes, outputs summary stats (to screen - saved in next function), creates and saves violin plots of
# nFeature_RNA, nCount_RNA and percent.mt to use in QC filtering, also creates and saves counts vs genes plot. Saves in Seurat directory in sample directory (must exist)
# Returns Seurat object. Flowcell argument for unique labelling of samples.
multisample_seurat_QC<-function(count_matrix, sample, flowcell, output1, output2) {
  if (length(grep(1, flowcell)) == 1) {
    z<-1
  } else {
    z<-2
  }
  
  unique_name<-paste(sample, "_FC", z, sep = "")
  
  seurat_obj<-CreateSeuratObject(counts = count_matrix, project = sample)
  
  seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  summary_stats <- data.frame(
    sample = sample,
    n_cells = ncol(seurat_obj),
    median_umis = median(seurat_obj$nCount_RNA),
    median_genes = median(seurat_obj$nFeature_RNA),
    median_percent_mt = median(seurat_obj$percent.mt, na.rm = T)
  )
  
  print(summary_stats)
  
  if (is.na(sum(seurat_obj$percent.mt))) {
    p1<- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3) + ggtitle(paste("QC for", sample))
  } else {
    p1<- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + ggtitle(paste("QC for", sample))
  }
  ggsave(paste(base, "/", flowcell, "/", sample, "/Seurat/", output1, sep = ""), width = 10, height = 6)
  
  p2<-FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(paste("Counts vs Genes for", sample))
  ggsave(paste(base, "/", flowcell, "/", sample, "/Seurat/", output2, sep = ""), width = 10, height = 6)
  
  return(seurat_obj)

}



### EXTRACT DATA ###
# Loop over all samples across both flowcells and create annotated count matrices, then create seurat objects from these (with samples differentiated by flowcell) and combine all in a list
seurat_list<-list()

for (x in 1:length(in_dirs)) {
    fc<-in_dirs[x]
    print(paste("Working on flowcell", x))
    
    for (i in 1:length(sample_names)) {
      sample<-sample_names[i]
      print(paste("Processing sample", sample, "number", i, "of", length(sample_names)))
      
      # Make sample specific output directory
      dir.create(paste(base, "/", fc, "/", sample, "/Seurat", sep = ""))
      
      # Run the functions made above
      counts<-make_annotated_matrix(fc, sample, gene_names)
      seurat_obj<-multisample_seurat_QC(counts, sample, fc, "seurat_prefilter_vln.png", "seurat_prefilter_count_v_gene.png")
      
      # Adjust for same samples on different flowcells - need unique names for list
      unique_name<-paste(sample, "_FC", x, sep = "")
      
      # Add flowcell column to seurat object so they can be differentiated later
      seurat_obj$flowcell<-x
      
      # Save sample specific seurat object
      saveRDS(seurat_obj, paste(base, "/", fc, "/", sample, "/Seurat/unprocessed_seurat.rds", sep = ""))
      
      # Add sample specif seurat object to list of all seurat objects
      seurat_list[[unique_name]]<-seurat_obj
    } 
  }

saveRDS(seurat_list, paste(base, "/outputs/seurat_processing/intermediate/seurat_list.rds", sep = ""))

# Create table of QC stats across all seurat objects in list
qc_stats <- lapply(names(seurat_list), function(sample) {
  data.frame(
    sample = sample,
    nCount_RNA = seurat_list[[unique_name]]$nCount_RNA,
    nFeature_RNA = seurat_list[[unique_name]]$nFeature_RNA,
    percent.mt = seurat_list[[unique_name]]$percent.mt
  )
}) %>% bind_rows()



### ANALYSE DATA ###

# View QC metric plots to identify suitable cutoffs for filtering 
count<-ggplot(qc_stats, aes(x=nCount_RNA, fill = sample)) + 
  geom_density(alpha=0.3) + 
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); count
ggsave(paste(base, "/outputs/seurat_processing/figures/qc/seurat_nCount_RNA.png", sep = ""), count, bg="white", width = 10, height = 6)

feature<-ggplot(qc_stats, aes(x=nFeature_RNA, fill = sample)) + 
  geom_density(alpha=0.3) +
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); feature
ggsave(paste(base, "/outputs/seurat_processing/figures/qc/seurat_nFeature_RNA.png", sep = ""), feature, bg="white", width = 10, height = 6)

mt<-ggplot(qc_stats, aes(x=percent.mt, fill = sample)) + 
  geom_density(alpha=0.3) +
  guides(shape = guide_legend(override.aes = list(size=0.5)), color = guide_legend(override.aes = list(size=0.5))) +
  theme(legend.title = element_text(size=3)) +
  theme(legend.text = element_text(size=3)); mt
ggsave(paste(base, "/outputs/seurat_processing/figures/qc/seurat_percent.mt.png", sep = ""), mt, bg="white", width = 10, height = 6)


# Set and apply filtering across all seurat objects in list
feat_cutoffs<-quantile(qc_stats$nFeature_RNA, probs=c(0.01, 0.99))
count_cutoffs<-quantile(qc_stats$nCount_RNA, probs=c(0.01,0.99))
mt_cutoff<-quantile(qc_stats$percent.mt, probs = 0.95, na.rm = TRUE)

if (mt_cutoff < 0.1){
  mt_cutoff <- 1
}

seurat_list_filtered<-seurat_list
for (unique_name in names(seurat_list)) {
  seu<-seurat_list[[unique_name]]
  
  seu<-subset(seu, 
              subset = 
                nFeature_RNA > feat_cutoffs[1] & 
                nFeature_RNA < feat_cutoffs[2] & 
                nCount_RNA > count_cutoffs[1] &
                nCount_RNA < count_cutoffs[2] &
                percent.mt < mt_cutoff)
  
  seurat_list_filtered[[unique_name]]<-seu
}

saveRDS(seurat_list_filtered, paste(base, "/outputs/seurat_processing/intermediate/seurat_list_filtered.rds", sep = ""))

# Identify how many cells were kept per sample
before<-sapply(seurat_list, function(x) ncol(x))
after<-sapply(seurat_list_filtered, function(x) ncol(x))

# Create comparison dataframe of n cells before and after filtering
cell_comparison<-data.frame(sample = names(seurat_list), before = before, after = after)

# Calculate number of cells removed by filtering for each sample
cell_comparison<- cell_comparison %>%
  mutate(removed = before - after)

# Add flowcell info
for (i in 1:length(cell_comparison$sample)) {
  cell_comparison$flowcell[i]<-tail(strsplit(cell_comparison$sample[i], "_")[[1]], n=1)
}

# Format dataframe ready to create stacked barplot
df_long<- cell_comparison %>%
  select(sample, flowcell, after, removed) %>%
  pivot_longer(cols = c("after", "removed"), names_to = "status", values_to = "count")

# Create stacked barplot of removed vs remaining cells per sample
# Save it
cell_status<-ggplot(df_long, aes(x=sample, y=count, fill=status))+
  geom_bar(stat="identity")+
  facet_wrap(~flowcell, scales = "free_x", ncol = 2)+
  labs(title="Number of Cells before and after filtering", x="Sample", y="Cell count")+
  theme_minimal()+
  scale_fill_manual(values = c("removed" = "darkred", "after" = "#33a02c"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)); cell_status
ggsave(paste(base, "/outputs/seurat_processing/figures/qc/seurat_filtering.png", sep = ""), cell_status, bg="white", width = 10, height = 6)

# Normalise each seurat_object in list
seurat_list_filtered_normalised<-lapply(seurat_list_filtered, function(x) {
  NormalizeData(x)
})

# Find variable features for each seurat object in list 
seurat_list_filtered_normalised_variable<-lapply(seurat_list_filtered_normalised, function(x) {
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# Save filtered, normalised, variable feature found list of seurat objects ready for integration and analysis
saveRDS(seurat_list_filtered_normalised_variable, paste(base, "/outputs/seurat_processing/intermediate/seurat_list_filtered_normalised_variable.rds", sep = ""))


