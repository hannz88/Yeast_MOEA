# R 3.6.1
# Optimization for PCA using genetic algorithm
# Performed using QMUL HPC Apocrita

# load dependencies----
# install library to local and getting the package from the repos
if(!require(genalg)){
  install.packages("genalg")
  library(rlist)
}
if(!require(FactoMineR)){
  install.packages("FactoMineR")
  library(FactoMineR)
}
if(!require(factoextra)){
  install.packages("factoextra")
  library(factoextra)
}

# load the file----
yeast_df <- read.csv("file path name", sep=",", header = TRUE, row.names = 1)
yeast_df <- na.omit(yeast_df)
yeast_df <- yeast_df[-5]  # omit the 5th column with cdc28


# the evaluation function for genalg----
evalfunc <- function(x){
  dataset <- yeast_df[x==0,]
  if (sum(x==1)>2000){
    return(42)
  }
  else {
    yeast_pca_ga <- PCA(dataset, scale.unit = TRUE, ncp = 2, graph = FALSE)
    eig_val_ga <- data.frame(get_eigenvalue(yeast_pca_ga))
    num_row <- nrow(eig_val_ga[eig_val_ga$cumulative.variance.percent<76,])
    return(eig_val_ga$cumulative.variance.percent[3]*-1)
  }
}

# design the model----
iter <- 500

ga_model <- rbga.bin(size = nrow(yeast_df), popSize = 200, iters = iter, 
                     elitism = TRUE, evalFunc = evalfunc)

# to get the best genome from the results----
best_genome <- vector()
best_fitness <- 0
df <- ga_model$population
for (i in 1:nrow(df)){
  x = df[i,]
  fitness <- evalfunc(x)
  if (fitness < best_fitness){
    best_fitness <- fitness
    best_genome <- x
  }
}

save.image("./yeast_3415_genalg.RData")

# explore the data
best_genome
new_df <- yeast_df[best_genome==0,]
nrow(new_df)
ncol(new_df)

new_pca <- PCA(new_df,scale.unit = TRUE, ncp=2, graph = TRUE)
fviz_pca_var(new_pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)   
fviz_eig(new_pca, addlabels = TRUE, ylim=c(0,45))
      