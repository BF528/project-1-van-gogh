library(gplots)

# Read in the data as a data frame
data <- data.frame(read.csv('/project/bf528/project_1/data/example_intensity_data.csv', sep=' ', header=TRUE))

############# 4. Noise Filtering & Dimensionality Reduction 
# 1) Expressed in at least 20% of samples (i.e. for each gene, at least 20% of the gene-expression values must be > log2(15)
expressed_filter <- data[rowSums(data > log2(15)) >= (0.2*ncol(data)), ]

# 2) Have a variance significantly different from the median variance of all probe sets using a threshold of p < 0.01
# Get the critical region
df<-ncol(expressed_filter)-1
chi_sq_lower<-qchisq(0.01/2,df)
chi_sq_upper<-qchisq(0.01/2,df,lower.tail = FALSE)
# Compute the test statistic for each gene: T = (N-1)(s/o)**2 where df=N-1, variance=s, median variance=o
expressed_filter$variance <- apply(expressed_filter, 1, var)
expressed_filter$test_stat <- df*expressed_filter$variance/median(expressed_filter$variance)
chi_filter<-subset(expressed_filter, test_stat > chi_sq_lower & test_stat < chi_sq_upper)

# 3) Have a coefficient of variation > 0.186
variation_filter<-subset(chi_filter, variance > 0.186)


# Write to a file
# Remove extra columns to write out to a file
variation_filter<-subset(variation_filter, select = -c(variance, test_stat))
chi_filter<-subset(chi_filter, select = -c(variance, test_stat))

# 4) Write out a different containing the gene expression matrix for genes passing all three filters
write.csv(variation_filter, '/projectnb/bf528/users/van-gogh/project_1/data/all_filtered_expression_matrix.csv')
# 5) Write out the expression matrix for probe sets that pass the expression threshold
write.csv(chi_filter, '/projectnb/bf528/users/van-gogh/project_1/data/probeset_filtered_expression_matrix.csv')

# The number of genes that pass all of these thresholds
nrow(variation_filter)

############### 5. Hierarchical Clustering & Subtype Discovery
# 1) Perform hierarchical clustering of patients
tree<-hclust(dist(t(variation_filter)))

# 2) Cut the dendogram so the samples are divided into 2 clusters
clusters<-cutree(tree, 2)
sum(clusters==1)
sum(clusters==2)

# 3) Create a heatmap of the gene-expression of each gene across all samples
# Load metadata
metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
metadata <- subset(metadata, select = c(geo_accession, cit.coloncancermolecularsubtype))
colors <- ifelse(metadata$cit.coloncancermolecularsubtype == 'C3', 'blue', 'red')

# Vector of colors for annotation: “red” if the sample belongs to the C3 subtype and “blue” otherwise
heatmap.2(as.matrix(variation_filter), xlab='Patient Sample', ylab='Gene', 
          main='Gene Expression of Genes Across Samples',
          ColSideColors = colors, trace='none', density.info = 'none',
          key.xlab='Gene Expression Level', margins=c(6,6))

# Legend for Molecular Subtype
legend(x=-0.05, y=0.9, xpd=TRUE, cex=0.75,
       legend=c('C3', 'Other'), title='Subtype', fill=c('blue', 'red'))