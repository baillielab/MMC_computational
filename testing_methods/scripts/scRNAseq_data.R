# Install the required packages
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("scRNAseq")

# Try to download/work with some data
library(scRNAseq)
test_data <- DarmanisBrainData(ensembl = FALSE, location = TRUE, remove.htseq = TRUE, legacy=FALSE)

# colData(test_data)$experiment_sample_name) is the donor