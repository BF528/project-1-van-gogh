# Project Description

A brief description of what this repository is for and what it contains

# Contributors
* Data Curator: Monil Gandhi 
* Programmer: Lindsay Wang (@LindsayW007)
* Analyst: Elysha Sameth (@esameth)
* Biologist: Andrew Gjelsteen (@agjelste)

# Repository Contents
## Programmer
### project1_PCA.R
### project1_norm_qc.R

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
