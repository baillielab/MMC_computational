rm(list=ls())

barcode_summary<-data.frame(sample=NA, FC=NA, barcode=NA, n_umis=NA, total_count=NA, rank=NA)

sample_names<-read.delim("/mnt/odap-beegfs/a008/pilot_pilot_new/data/samples.tsv")
sample_names<-sample_names[,1]
if (length(grep("Saccular", sample_names)) != 0) {
  sample_names<-sample_names[-c(grep("Saccular", sample_names))] # Removing saccular as not part of MMC analysis
}

in_dirs<-c("/mnt/odap-beegfs/a008/pilot_pilot_new/outputs/arcas_pipeline/FC1_unpaired/MMC", 
           "/mnt/odap-beegfs/a008/pilot_pilot_new/outputs/arcas_pipeline/FC2_unpaired/MMC")

knee_threshold<-100

for (j in 1:length(in_dirs)) {
  # Run through each input dir/flowcell
  fc<-in_dirs[j]
  
  for (i in 1:length(sample_names)) {
    # Process each sample
    sample<-sample_names[i]
    
    barcode_sample<-read.delim(paste(fc, "/", sample, "/QC/barcode_summary.tsv", sep = ""))
    
    barcode_sample$sample<-sample
    barcode_sample$FC<-j
    
    barcode_sample<-barcode_sample %>%
      arrange(desc(n_umis)) %>%
      mutate(rank = row_number())
    
    barcode_summary<-rbind(barcode_summary, barcode_sample)
    
    filtered_barcodes <- barcode_sample %>%
      filter(n_umis >= knee_threshold) %>%
      pull(barcode)
    
    write.csv(filtered_barcodes, paste(fc, "/", sample, "/QC/barcode_filter.csv", sep = ""), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  }
}

FC2<-barcode_summary[which(barcode_summary$FC == 2),]
FC1<-barcode_summary[which(barcode_summary$FC == 1),]

    
# Generate knee plot with knee rank annotated
# Save plot
knee1<-ggplot(FC1, aes(x=rank, y=n_umis, colour = sample)) +
  geom_line() +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Barcode Rank Plot (log-log)", x = "Rank", y = "UMI count")

knee1

knee2<-ggplot(FC1, aes(x=rank, y=n_umis, colour = sample)) +
  geom_line() +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Barcode Rank Plot (log-log)", x = "Rank", y = "UMI count")

knee2

