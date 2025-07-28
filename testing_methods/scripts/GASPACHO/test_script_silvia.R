#!/usr/bin/env Rscript

# The data used in this script is provided by the GASPACHO github repo - https://github.com/natsuhiko/GASPACHO/tree/v1.0.0
# There are 68 individuals with number of cells ranging from 5 to 674 cells. 

metadata = readRDS("../data/GASPACHO/metadata.RDS")
cpm = readRDS("../data/GASPACHO/log_cpm_4999_22188.RDS") 
init_param = readRDS("data/init_param.RDS") # parameters for GASPACHO

# Inspect the files and understand the structure

str(metadata) 

# The metadata is a data.frame

gplvm = GPLVM(cpm, metadata,
	Xi     = init_param$Xi, # latent variables
	Ta     = init_param$Ta, # inducing variables
	delta0 = init_param$delta, # variance parameters for metadata
	Alpha0 = init_param$Alpha, # fixed effect (NULL)
	sigma0 = init_param$sigma2, # residual variance for genes
	omega0 = init_param$omega2, # residual variance for cells
	zeta0  = init_param$zeta, # grand mean
	rho0   = init_param$rho, # length parameters for Kernels
	theta0 = init_param$theta, # variance parameter for GP
	ITRMAX = 1000)

# Estimating the Donor by Context interaction effect
gplvm = updateDxCSE(as.matrix(cpm), gplvm)

# Computing eQTL Bayes factors
G_OAS1 = readRDS("Data/G_OAS1.RDS")
# yj : expression vector for a target gene
# G  : genotype dosage matrix (locus x donor). Non-normalised (each value in [0,2])
# did: Donor ID compatible with the column of G
bfs = getBF(yj = as.numeric(cpm[3329,]), gplvm, G = as.matrix(G_OAS1[,10:77]), 
            did = as.numeric(as.factor(metadata$donor)))