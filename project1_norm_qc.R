#!/usr/bin/env Rscript

#setwd("//projectnb/bf528/users/van-gogh/project_1")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install(c("affy", "affyPLM", "sva", "AnnotationDbi", "hgu133plus2.db"))

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

#list the name of the files
data_list <- list.celfiles("/projectnb/bf528/users/van-gogh/project_1/samples")
#Extract sample names from the file names
sample_name <- sapply(strsplit(data_list, "_"), '[', 1)
#data_list <- dir(pattern = "CEL.gz")
array <- ReadAffy(filenames=data_list, verbose = TRUE, sampleNames=sample_name, celfile.path="/projectnb/bf528/users/van-gogh/project_1/samples")

#Normalization and Quality Control
eset <- rma(array, subset=NULL, normalize=TRUE, background=TRUE)
#Convert AffyBatch into PLM object
Pset <- fitPLM(array, normalize=TRUE, background=TRUE)
#Calclate RLE statistics
rle_stat <- RLE(Pset, type="stats")
#Calclate NUSE statistics
nuse_stat <- NUSE(Pset, type="stats")
#Histogram for RLE median
hist(rle_stat[1,], xlab="median RLE scores", main="Historgam of Median RLE scores")
#Historgam for NUSE median
hist(nuse_stat[1,], xlab="median NULE scores", main="Historgam of Median NULE scores")

##Batch Correction
annotation_file <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv", header=TRUE)
combat_var <- annotation_file[c("geo_accession", "normalizationcombatbatch", "normalizationcombatmod")]
#Set the index to be the $geo_accession
rownames(combat_var) <- combat_var$geo_accession

pheno <- pData(eset)
batch <- annotation_file$normalizationcombatbatch
mod <- model.matrix(~as.factor(annotation_file$normalizationcombatmod), data=pheno)
norm_data <- ComBat(dat=exprs(eset), batch=batch, mod=mod)

write.csv(norm_data, "/projectnb/bf528/users/van-gogh/project_1/norm_data.csv")
