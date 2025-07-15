# Single-cell eQTL mapping - project specification

**Author(s):** Silvia  
**Last updated:** 14-07-2025  
**Other contributors:**  

## 1. Project overview

**Aims/hypotheses**

Aim 1: Which genetic variants are associated with differences in gene expression levels in the population?  
Hypothesis tested: SNP *x* is statistically associated with expression of gene *Y*.

    A) Pseudo-bulk method (inspired by 1k1k):  @@@Check where the covariates are pre-corrected for, either after step 1 or step 2
    - Cluster into cell types - static, using PCA or other methods.  
    - Calculate gene expression values (one count per individual per cell type).  
    - Precorrect gene expression values for technical and population covars (I *think*)
    - Calculate correlation (Spearman's rank) between SNP and residuals from precorrection to identify lead SNPs associated with changes in expression.  
    - Correct gene expression levels using linear model for SNP.  
    - Calculate correlation between residuals and remaining SNPs (do this and previous step up to 5 times to get conditional associations). 

    B) 'Dynamic methods' @@@More research required
    - Model cell states dynamically
    - Mixed-model association 

**Key outputs**

A list of association summary statistics for each gene, eg:

|Gene | eQTLs (SNPs) | p-values | other stats|
|----|----|----|----|
| XYZ | rsxxxx or chr:loc | 5 x 10^-12 | beta? spearman's rank coeff (rho?)|

## 2. Functional & design specification

**Input data**

A) Pseudo-bulk method:
- Count matrix
- Covariates for each 

**intermediary data**

[gene x individual] for each cell

**Pre-processing steps**

**eQTL design choices**
- What method is used to produce cell-types?
- What method/model is used to assess association?
- What covariates are used?
- Are we doing conditional analysis or joint analysis?

## Analysis plan and timeline

| Task                        | Week 1 (Jul 14) | Week 2 (Jul 21) | Week 3 (Jul 28) | Week 4 (Aug 4) | Week 5 (Aug 11) |
|-----------------------------|------------------|----------------|-----------------|------------------|------------------|
| Pipeline design             | █████████████████|████████████████|█████████████████|                  |                  |
| Pipeline building           |                  |████████████████|█████████████████|██████████████████|                  |
| Test on ODAP                |                  |                |                 |                  |                  |
| Test on parse/1k1k data     |                  |                |█████████████████|                  |                  |

## Open questions
