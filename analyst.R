library(gplots)

# Read in the data as a data frame
data<-data.frame(read.csv('/project/bf528/project_1/data/example_intensity_data.csv', sep=' ', header=TRUE))

#################################################
# 4. Noise Filtering & Dimensionality Reduction 
#################################################
# 1) Expressed in at least 20% of samples (i.e. for each gene, at least 20% of the gene-expression values must be > log2(15)
expressed_filter<-data[rowSums(data > log2(15)) >= (0.2*ncol(data)), ]

# 2) Have a variance significantly different from the median variance of all probe sets using a threshold of p < 0.01
# Perform an upper one-tailed alternative
df<-ncol(expressed_filter)-1
chi_sq_upper<-qchisq(1-0.01, df, lower.tail=FALSE)
# Compute the test statistic for each gene: T = (N-1)(s/o)**2 where df=N-1, variance=s, median variance=o
expressed_filter$variance <- apply(expressed_filter, 1, var)
expressed_filter$test_stat <- df*expressed_filter$variance/median(expressed_filter$variance)
chi_filter<-subset(expressed_filter, test_stat > chi_sq_upper)
# Remove extra columns to write out to a file
chi_filter<-subset(chi_filter, select = -c(variance, test_stat))

# 3) Have a coefficient of variation > 0.186
variation_filter<-subset(chi_filter, apply(chi_filter, 1, function(x) sd(x)/mean(x)) > 0.186)


#################################################
# Write to a file
#################################################
# 4) Write out a different containing the gene expression matrix for genes passing all three filters
write.csv(variation_filter, '/projectnb/bf528/users/van-gogh/project_1/data/all_filtered_expression_matrix.csv')
# 5) Write out the expression matrix for probe sets that pass the expression threshold
write.csv(chi_filter, '/projectnb/bf528/users/van-gogh/project_1/data/differential_filtered_expression_matrix.csv')

# The number of genes that pass all of these thresholds
print(paste0('Number of genes that pass all thresholds: ', nrow(variation_filter)))


#################################################
# 5. Hierarchical Clustering & Subtype Discovery
#################################################
# Function to perform the clustering because it needs to be done twice: for filter 2 and after filter 3
welch_t_test<-function(data_frame, file_name) {
  # 1) Perform hierarchical clustering of patients
  tree<-hclust(dist(t(data_frame)))
  
  # 2) Cut the dendogram so the samples are divided into 2 clusters
  clusters<-cutree(tree, 2)
  print(paste0('Number of samples in cluster 1: ', sum(clusters==1)))
  print(paste0('Number of samples in cluster 2: ', sum(clusters==2)))
  
  # 4) identify genes differentially expressed between the two clusters using a Welch t-test
  cluster1<-data_frame[, clusters==1]
  cluster2<-data_frame[, clusters==2]
  
  # Create an empty list to fill with dataframes
  t_test<-vector('list', nrow(data_frame))
  # For every gene...
  for (row in 1:nrow(data_frame)) {
    # Perform a Welch t-test for cluster1 and cluster2
    welch<-t.test(cluster1[row,], cluster2[row,])
    # Create a list of dataframes with the columns 't', 'p' 
    t_test[[row]]<-data.frame(t=welch$statistic, p=welch$p.value)
  }
  
  # Bind the list of dataframes into one dataframe
  t_test<-do.call(rbind, t_test)
  # Set the gene names as the row names
  row.names(t_test)<-row.names(data_frame)
  # Get padjust
  t_test$padj<-p.adjust(t_test$p, method='fdr')
  # Write out to a comma separated file 
  write.csv(t_test, file_name)
  
  # Print statements
  padj_0.05<-subset(t_test, padj < 0.05)
  padj_0.05<-padj_0.05[order(padj_0.05$padj),]
  print(paste0('Number of genes differentially expressed at padj < 0.05 between the clusters: ', nrow(padj_0.05)))
  print('The most differentially expressed genes that best defines each cluster were determined by the criteria supplied by the supplementary: adjusted p-value < 1e-5 and |log2 fold change| > 0.5')
  print('The top differentially expressed genes based on adjusted p-value < 1e-5 are: ')
  print(paste0(head(rownames(subset(padj_0.05, padj < 1e-5)))))
}

# Perform the t-test analysis
print('T-test for expression matrix with all filters')
welch_t_test(variation_filter, '/projectnb/bf528/users/van-gogh/project_1/data/all_filters_t_test.csv')
print('T-test for expression matrix for probesets that pass the expression threshold')
welch_t_test(chi_filter, '/projectnb/bf528/users/van-gogh/project_1/data/differential_filter_t_test.csv')

# 3) Create a heatmap of the gene-expression of each gene across all samples
# Load metadata
metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
metadata <- subset(metadata, select = c(geo_accession, cit.coloncancermolecularsubtype))
metadata <- metadata[metadata$geo_accession %in% colnames(variation_filter),]
colors <- ifelse(metadata$cit.coloncancermolecularsubtype == 'C3', 'blue', 'red')

# Save heatmap as a file
png('/projectnb/bf528/users/van-gogh/project_1/data/heatmap.png', width=1920, height=1080, res=100)
# Vector of colors for annotation: “red” if the sample belongs to the C3 subtype and “blue” otherwise
heatmap.2(as.matrix(variation_filter), xlab='Patient Sample', ylab='Gene', 
          main='Gene Expression Across Samples',
          ColSideColors = colors, trace='none', density.info = 'none',
          key.xlab='Expression Level', scale='row', margins=c(7,7))

# Legend for Molecular Subtype
legend(x=0.9, y=1, xpd=TRUE, inset=c(-0.15,0),
       legend=c('C3', 'Other'), title='Molecular Subtype', fill=c('blue', 'red'))
dev.off()
