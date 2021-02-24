setwd('/projectnb/bf528/users/van-gogh/project_1/data')
data <- read.csv('differential_filter_t_test.csv', header=TRUE)
# install.packages("xlsx")
library(xlsx)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")
# BiocManager::install(pkgs = c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("reshape")

library(GSEABase)
library(reshape)
library(dplyr)
library(tidyverse)
library(hgu133plus2.db)
library(RSQLite)
library(AnnotationDbi)


data <- rename(data, c("PROBEID" = "X"))

matches <- AnnotationDbi::select(hgu133plus2.db, 
                                 keys = as.character(data$PROBEID), 
                                 columns = ("SYMBOL"))
#returns gene symbols for the corresponding probeIDs (Column 1 of the .csv file)

# Remove duplicates based on $PROBEID column name, .keep_all is TRUE to not remove it from data set:
res.list <- matches
curatedMatches <- res.list %>% distinct(PROBEID, .keep_all = TRUE)

newData <- merge(data, curatedMatches, by = "PROBEID")


##### End of combing through data ######
## select the top 1000 up- and 1000 down-regulated (i.e. positive and negative -log2 fold change) genes, irrespective of significance.
## generating table with the top 10 up- and down- regulated genes ##


sortds1_up <- newData[newData$t > 0,] #subset with only + t values.
sortds2_up <- sortds1_up[order(sortds1_up$padj), ] #sorts by ascending padj
sortds3_up <- sortds2_up[1:1000,]
sortds4_up <- sortds3_up[order(sortds3_up$padj), ]
sortds5_up <- sortds4_up[1:10,] #Selecting top 10 most up-regulated

sortds1_down <- newData[newData$t < 0,] #Subset with only negative t values
sortds2_down <- sortds1_down[order(sortds1_down$padj),] #sort by ascending padj
sortds3_down <- sortds2_down[1:1000,]
sortds4_down <- sortds3_down[order(sortds3_down$padj), ]
sortds5_down <- sortds4_down[1:10,]#Selecting top 10 most down-regulated

sortTable <- rbind(sortds5_up, sortds5_down)
sortTable
# write.xlsx(sortTable, "/projectnb/bf528/users/van-gogh/project_1/differential_expr_table.xlsx") #for exporting it

##### Read in .gmt files and use them for analysis: #####

c2gmt <- getGmt("./c2.cp.kegg.v7.2.symbols.gmt",collectionType=BroadCollection(category="c2"),geneIdType=SymbolIdentifier()) #186 elements
c5gmt <- getGmt("./c5.go.v7.2.symbols.gmt",collectionType=BroadCollection(category="c5"),geneIdType=SymbolIdentifier()) #186 elements
hgmt <- getGmt("./h.all.v7.2.symbols.gmt",collectionType=BroadCollection(category="h"),geneIdType=SymbolIdentifier()) # 10271 elements

# Write a function that accepts a set of differentially expressed gene symbols and a single gene set to test code
# Return contingency table to pass to the Fisher.test function

sortds4_up #the up-regulated to be used here
sortds4_down #the down-regulated to be used here:


geneIds <- geneIds(c2gmt)
geneIds
hi <- (geneIds[1])
# sapply(hi[1], length)
# y <- 0 #the index to get the name of a gene set (up-regulated)
# z <- 0 #the index to get the name of a gene set (down-regulated)
# name_up <- list()
# name_down <- list()

#### Creating contingency table: ####
# nrow(newData)





common_Genes <- function(geneCollection, geneSet) {   #takes in geneCollection (.gmt) and geneSet
  y <- 0
  geneIds <- geneIds(geneCollection)
  for (i in geneIds) {
    table_list <- list
    y <- y + 1 #integer to serve as token to find name of Gene Set using following line:
    name <- names(geneIds[y]) #The name of each individual Gene Set from geneCollection
    
    
    common <- intersect(geneSet$SYMBOL,i)
    commonDiffGenes <- length(intersect(geneSet$SYMBOL,i)) #1,1 in table
    uncommonGenes <- nrow(geneSet) - (commonDiffGenes) #1, 2 in table
    totalDiffExpr <- (uncommonGenes + commonDiffGenes) #1,3 in table
    
    geneSetCommonDiff <- length(intersect(newData$SYMBOL,i)) - commonDiffGenes #(2,1) in table
    geneSetUncommonDiff <- nrow(newData) - geneSetCommonDiff #(2,2)
    notDiffTotal <- geneSetCommonDiff + geneSetUncommonDiff #(2,3)
    
    totalGeneInSet <- commonDiffGenes + geneSetCommonDiff #(3,1)
    totalOutOfSet <-  uncommonGenes + geneSetUncommonDiff #(3,2)
    totalGenes <- totalGeneInSet + totalOutOfSet #(3,3)
    
    fisher_result <- fisher.test(matrix(c(commonDiffGenes,uncommonGenes,geneSetCommonDiff,geneSetUncommonDiff),nrow=2))
    
    table_GeneCollection <- data.frame(x1 = 0,                      # Create example data
                                       x2 = 0)
    new <- rep(i, ncol(data))                       # Create new row
    data[nrow(data) + 1, ] <- new                   # Append new row
    
    
    if (length(common) > 10){
      ### create contingency table:
      A <- c(commonDiffGenes, uncommonGenes, totalDiffExpr)
      B <- c(geneSetCommonDiff, geneSetUncommonDiff, notDiffTotal)
      C <- c(totalGeneInSet, totalOutOfSet, totalGenes)
      contingencyT <- data.frame("differentially-expressed" = A, "not differentially-expressed" = B, "Total" = C)
      # return(contingencyT)
    }
    return (fisher_result$p.value)
  }
  return ()
}

sortds4_down$fisher.p.value <- common_Genes(c2gmt, sortds4_down)


# p.adjust(method = "FDR")
