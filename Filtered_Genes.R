setwd("/Users/bigyambat/Desktop/RNASeq_Project")

library(readxl)
library(ggfortify)

Filtered_Genes = read_excel('Bigy-gene_exp.xlsx', 'Sheet1')

#Inner Join
Merged_Filtered =  merge(x=Filtered_Genes, y=Master_FPKM, by= "gene_id")

Merged_Filtered2 = Merged_Filtered[, -c(1:14)]


#Removing Sum

#Merged_Filterd_No0Sum = subset(Merged_Filtered2, Sum!=0)

Merged_Filter_rmSum <- Merged_Filtered2[1:(length(Merged_Filtered2)-1)]

Merged_Filter.rmSum2 <- Merged_Filter_rmSum[-c(4087), ]

#Creating Pseudo count of 0.1 and log(10) conversion
Merged_Filter.psuedo <- Merged_Filter.rmSum2[7:14] + 0.1

Merged_Filter.log <- log10(Merged_Filter.psuedo) #log change

Merged_Filter.log_naomit <- na.omit(Merged_Filter.log)

Merged_Filter.log_naomit2 <- Merged_Filter.log_naomit[-c(4087), ]

Merged_Filter_Final <- cbind(Merged_Filter.rmSum2[1:7], Merged_Filter.log_naomit2) 


#PCA Time

transposed_Filter <- t(Merged_Filter.log_naomit2)

transposed_Filter <- as.data.frame(transposed_Filter)

transposed_Filter$Sample <- rownames(transposed_Filter)

transposed.data <- transposed_Filter[c(1:190)]

Filter_PCA <- autoplot(prcomp(transposed.data), data = transposed_Filter, colour = 'Sample', )

print(Filter_PCA + ggtitle("PCA Plot using Filtered Genes"))


#Dendrogram

scaled.transposed.data <- scale(transposed_Filter[c(1:190)])

head(scaled.transposed.data, n=3)

dist.eucl <- dist(scaled.transposed.data, method = "euclidean")

#dist.eucl$Samples <- row.names(scaled.transposed.data)

library(factoextra)

dist.cor <- get_dist(scaled.transposed.data, method = "pearson")


hclust.ward.eucl <- hclust(d = dist.eucl, method = "ward.D2")
hclust.ward.cor <- hclust(d = dist.cor, method = "ward.D2")
hclust.complete.eucl <- hclust(d = dist.eucl, method ="complete")
hclust.complete.cor <- hclust(d = dist.cor, method ="complete")

plot(hclust.ward.eucl, labels =FALSE, main ="Euclidian - Ward's")

plot(hclust.ward.cor, labels = FALSE, main ="Correlation - Ward's")

plot(hclust.complete.eucl, labels = FALSE, main ="Euclidian - Complete")

plot(hclust.complete.cor, labels = FALSE, main ="Correlation - Complete")

#Heatmaps

library(tidyverse)

df_test <- scale(mtcars)

Heatmap1 <- Merged_Filter_Final[-c(1:3, 5,6)]

Heatmap2 <- Heatmap1[!duplicated(Heatmap1$gene_short_name), ]

rownames(Heatmap2) <- Heatmap2[,1]

Heatmap3 <-Heatmap2[-c(1, 2)]

Heatmap4 <-data.matrix(Heatmap3)

heatmap(Heatmap4, scale = "none")




