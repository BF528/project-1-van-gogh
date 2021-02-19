#!/usr/bin/env Rscript

library(ggplot2)

ds <- read.csv("norm_data.csv", row.names=1, header=TRUE)
ds_pca <- t(scale(t(ds)))
pca <- prcomp(ds_pca, scale=FALSE, center=FALSE)

metadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv", header=TRUE)

#plot PCA
pca_df <- data.frame(PC1=pca$rotation[,1], PC2=pca$rotation[,2], Cancer.Subtype=metadata$SixSubtypesClassification)
#png("PCA_Result.png")
ggplot(pca_df, aes(PC1, PC2, col=Cancer.Subtype, fill=Cancer.Subtype)) + 
  stat_ellipse(geom='polygon', col='black', size=0.2, alpha=0.2) + 
  geom_point(shape=21, col='black') + 
  xlab(paste("PC1 -", summary(pca)$importance[2,1]*100, "%")) + 
  ylab(paste("PC2 -", summary(pca)$importance[2,2]*100, "%")) +
  labs(fill="Subtypes") + 
  theme_bw() + ggtitle("PCA Result")
#dev.off()
