# Project Description

This analysis focuses only on reproducing the results from the comparison of the C3 and C4 tumor subtypes from Marisa, et al. The study was conducted in a two-phase design, where an initial set of “discovery” samples was used to identify patterns among the samples, and a separate set of “validation” samples was used to test if the results from the discovery set were robust. For this analysis, we have combined the discover and validation set samples into a single dataset. There are 134 samples in total.  
In the reproduction of the Marisa, et al. study, we retrieved gene expression values from GEO, normalized and conducted gene quality analysis, performed hierarchical clustering of samples, and determined significantly enriched gene sets.

# Contributors
* Data Curator: Monil Gandhi (@gandhimonil9823)
* Programmer: Lindsay Wang (@LindsayW007)
* Analyst: Elysha Sameth (@esameth)
* Biologist: Andrew Gjelsteen (@agjelste)

# Repository Contents
## Analyst
### analyst.R
To run this, use the command R analyst.R

* Reads in `norm_data.csv` created by `project1_norm_qc.R`
* Filters genes based on three metrics defined by Marisa, et al.
* Writes a gene expression matrix for genes passing all three filters to `all_filtered_expression_matrix.csv`
* Write a gene expression matrix for genes with a variance significantly different from the median variance of all probe sets to `differential_filtered_expression_matrix.csv`
* Performs unsupervised hierarchical clustering of C3 and C4 tumor subtypes
* Writes the results of Welch's t-test between the clusters and their p-values and FDR adjusted p-values to `all_filters_t_test.csv` (all three filters) and `differential_filter_t_test.csv` (filter two only)
* Outputs a heatmap 

## Biologist
### biologistScript.R
To run this, use the command R biologistScript.R while in the same directory as the differential_filter_t_test.csv as well as KEGG, GO and Hallmark gene set .gmt files.

* Reads in `differential_filter_t_test.csv` created by `analyst.R`
* Matches gene symbols to corresponding probe IDs using hgu133plus2.db and adds symbols to `differential_filter_t_test.csv` as a new column
