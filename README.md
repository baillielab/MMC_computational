# Single-cell molecular QTL mapping

Welcome to the downstream analysis arm of the Molecular Mechanisms Cluster (MMC) at the University of Edinburgh. If you have any questions, ideas or would like to collaborate, please reach out to silvia.shen[at]ed.ac.uk.  

**Aims**

Our aims are as follows (taken from the [FGI website](https://www.ukfunctionalgenomics.com/research/molecular-mechanisms-cluster/)). 

"This research cluster will tackle the “missing link” between the genome and disease. Our aim is to identify molecules in our bodies that cause disease by looking at how tiny variations in the genome change these molecules in single cells. We will create a unique dataset by reading molecular signals in human cells donated by hundreds of patients, potentially leading to ideas for new, effective drugs. ... In the Molecular Mechanisms Cluster we hold ourselves to a simple measure of success: the number of disease genes explained (number of disease-associated common genetic variants significantly colocalising with molecular quantitative trait loci)."

In simple terms, our goals are:  
⭐️ Identify QTLs  
⭐️ Create open source database of QTLs  
⭐️ Colocalise QTLs  
⭐️ Further analysis of QTLs  

For an overview of the whole process, please see the [MMC powerpoint](https://uoe-my.sharepoint.com/:p:/r/personal/kcampb2_ed_ac_uk/Documents/MMC_PPA_notes.pptx?d=w3afaa31d5ccd4efea7dff2bf2b120465&csf=1&web=1&e=jjdB1P). For an overview of what we could do with the very exciting single-cell data, see figure below (the projects I am currently working on are highlighted in red).   

**Key resources**  

- To see how key analysis decisions were made, please refer to the [specification documents](specification_docs).  
`<link to key resources that we can share>`

## Cell type analysis

Please see the [cell-type analysis plan](specification_docs/cell_type_plan.md) for what we plan to do. 

- What known cell types are there in our samples?
- What cells are communicating with each other? (receptor-ligand pairs)
- How does the abundance of cell-types change with environmental perturbations?

## eQTL mapping  
  
Please see the [single-cell eQTL mapping strategy](specification_docs/sceQTL_plan.md) on proposed options for single-cell eQTL mapping.  
For pipeline design, please see the [pipeline specification doc](pipeline_specification.md).  

## Future analysis ideas  
  
Please feel free to add any and all ideas you would like to implement (or see implemented by the team/me)!  
  
- Variance QTL mapping  
- Dynamic cell state/continuous cell state or phenotype mapping (perhaps go straight into this?)  
- Allelic imbalance analysis (see [key literature](specification_docs/key_literature.md))  
- Something with reads or tails of reads?   
- Producing a sc-eQTL atlas (similar to cattle gtex)  
- Which genes have the most variation across a single cell type? (and which SNPs/regions is this associated with) --> varQTLs
- Single-cell integration with polygenic risk scores? [Paper](https://www.nature.com/articles/s41588-022-01167-z)
  
![Image](specification_docs/potential_research_avenues.png)

## Resources

### General single-cell stuff
- The Theis lab's [book](https://www.sc-best-practices.org/preamble.html) on single-cell best practices and the accompanying [paper](https://www.nature.com/articles/s41576-023-00586-w) and [tutorial](https://www.embopress.org/doi/full/10.15252/msb.20188746) and [github](https://github.com/theislab/single-cell-tutorial).
- Aarun Lun's thoughts on single cell analysis on [github](https://ltla.github.io/SingleCellThoughts/)
- Sanger's Hemberg group resources: [book](https://www.singlecellcourse.org/)
- [scanpy tutorials](https://scanpy.readthedocs.io/en/latest/index.html)

### Single-cell QTL mapping resources

- Cuomo's [review paper: Single-cell genomics meets human genetics](https://www.nature.com/articles/s41576-023-00599-5). Great initial overview of main methods and challenges.
- Maria's [review paper: The Power of Single-Cell RNA Sequencing in eQTL Discovery](https://www.mdpi.com/2073-4425/13/3/502). Great starting point, v high level. 
- Zhang and Zhao's [review paper: eQTL studies: from bulk tissues to single cells](https://www.sciencedirect.com/science/article/pii/S1673852723001133). Has more equations and is more detailed on the QTL mapping bit.
