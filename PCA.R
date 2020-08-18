# R 3.6.1
# PCA analysis of the data with 3415 rows
# 24 July, 2020

# load dataset----
data <- read.csv("~/file_path_name", 
                 sep=",", header = TRUE, row.names = 1)
data <- na.omit(data)
data <- data[-5]  # omit the fifth column with cdc28
data <- scale(data)  # standardize

# load dependencies----
library(FactoMineR)
library(factoextra)

# compute PCA
pca_yeast <- PCA(data, scale.unit = TRUE, ncp=2, graph = TRUE)

# scree plot
fviz_eig(pca_yeast, addlabels = TRUE, ylim=c(0,30))

# positively correlated are grouped together
# variables far away are well-represented by the PCA  
fviz_pca_var(pca_yeast, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) 

# eigenvalue
pca_yeast$eig