rm(list=ls())

args<-commandArgs(trailingOnly = TRUE)

### IMPORT THINGS ###
# Load libraries
.libPaths("/mnt/odap-beegfs/software/R-vm/R/x86_64-pc-linux-gnu-library/4.5")

library(BUSpaRse)
library(Biostrings)

# Set file paths
base=args[1]

gtf<-paste(base, "/data/reference/homo_sapiens/GRCH38/reference.cdna.gtf.gz", sep = "")
fasta<-readDNAStringSet(paste(base, "/data/reference/homo_sapiens/GRCH38/reference.dna.fa.gz", sep = ""))

names(fasta)<-sub(" .*", "", names(fasta))


t2g<-tr2g_gtf(file = gtf,
              Genome = fasta)


t2g$transcript<-sub("\\..*", "", t2g$transcript)
t2g$gene<-sub("\\..*", "", t2g$gene)


write.csv(t2g, paste(base, "/data/reference/t2g.csv", sep = ""), row.names = FALSE, quote = FALSE)

