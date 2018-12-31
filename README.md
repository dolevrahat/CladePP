# CladePP

**Introduction**

The code enclosed in this repository provides an extension of the Normalized Phylogenetic Profiling algorithm describe in Tabach, Yuval, et al. "Identification of small RNA pathway genes using patterns of phylogenetic conservation and divergence." Nature 493.7434 (2013): 694.

Phylogenetic profiling (PP) is an evolution base method for prediction of gene function or interaction. PP relies on the assumption that genes that are functionally related (for example, participate in the same cellular pathway or are members of a common protein complex) manifest correlated patterns of occurrence and absence across the phylogenetic tree. 
In PP analysis, each gene is represneted using a phylogeneitc profile, a vector indiacating the presence or absence of an homolog for the gene across many evoluniarily diverse species.
Given a gene or a set of genes of interest, the algorithm searches for additonal genes with similiar phylogenetic profiles, thus highlighting potential funcational or interaction partners of the query gene(s).

While most current PP methods examine co-evolution between genes across the entire tree of life (global co-evolution), CladePP combines signals from global co-evolution with signals from gene pairs that are co-evolved within the context of a specific clade (local co-evolution) in order to refine predictions.

**Perquisites**

CladePP is implemented in R and requires the R libriries plyr, dplyr and optparse avilable from CRAN. Calculating a false discovery rate for the predictions also requires the qvalue library from Bioconductor.

**Generating CladePP predictions**


