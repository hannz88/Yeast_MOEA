# R 3.6.1
# post processing of the NSGA2.R output

# load dependencies----
if(!require(rlist)){
  install.packages("rlist")
  library(rlist)
}
if(!require(ecr)){
  install.packages("ecr")
  library(ecr)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(plotly)){
  install.packages("plotly")
  library(plotly)
}
if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}
if(!require(locStra)){
  install.packages("locStra")
  library(locStra)
}
if(!require(dict)){
  install.packages("dict")
  library(dict)
}
require(yeastfitness)

# load the file that has the functions to use nsga2 as white box----
source("file path for NSGA2_subfunctions.R")

# load RData----
load("file path for latest RData output from NSGA.R")


# obtain the Pareto optimal set----
fit_test <- as.matrix(pool$fitness)  # get the objectives of each individual
ranks <- doNondominatedSorting(fit_test)  # perform a non-dominated sorting on them
pareto_front_genome <- pool$population[ranks$ranks==1]
length(pareto_front_genome)  # this is Pareto optimal 

# Generate heatmaps for all of them----
for (i in 1:length(pareto_front_genome)) {
  target_genome <- pareto_front_genome[[i]]
  row_vect <- target_genome[1:nrow(data)]
  col_vect <- target_genome[-(1:nrow(data))]
  bic <- data[c(row_vect==1), c(col_vect==1)]
  paths <- paste0("file path/heatmap",i, ".png")
  png(filename = paths, width=600)
  heatmap.2(as.matrix(bic), col = redgreen(75), scale="row", key=T, keysize=1.5, density.info="none",trace="none",cexCol=0.9)
  dev.off()
}

# using jaccard coeeficients to compare the similarity and combine the solutions----
bin_list <- pool$population
# make a jaccard matrix
# turn the genome into a matrix with each genome occupying one column
genome_matrix <- matrix(unlist(bin_list), ncol=length(bin_list), byrow=FALSE)
sparse_genome_matrix <- Matrix(genome_matrix, sparse = TRUE)
# get a jaccard matrix between all the genome
jac_genome_matrix <- jaccardMatrix(sparse_genome_matrix)

# combine the genome with > 0.95 similarity
create_new_genome <- function(genome, sim_matrix){
  new_genome <- list()
  joined <- numvecdict()
  found = 0
  for (i in 1:nrow(sim_matrix)) {
    gene <- sim_matrix[i,]
    gene[gene==1.0] <- 0.0
    max_similarity <- max(gene)
      if (max_similarity >= 0.95) {
      found <- found + 1
      max_index <- which(gene==max_similarity)
      if ((i %in% joined$keys()) & (max_index %in% joined[[i]])) {
        next
      }
      combine <- as.integer(genome[[i]] + genome[[max_index]] > 0)
      new_genome <- list.append(new_genome, combine)
      if (max_index %in% joined$keys()) {
        joined$append_number(max_index, i)
      } else {
        joined$append_number(max_index, i)
      }
    } else {
      new_genome <- list.append(new_genome, genome[[i]])
    }
  }
  results <- list(combined_genome = new_genome, total_solutions = found)
  return(results)
}

tmp <- create_new_genome(genome = bin_list, sim_matrix = jac_genome_matrix)
found <- tmp$total_solutions
while (found) {
  print(length(tmp$combined_genome))
  new_bin_list <- tmp$combined_genome
  new_genome_matrix <- matrix(unlist(new_bin_list), ncol=length(new_bin_list), byrow=FALSE)
  new_sparse_genome_matrix <- Matrix(new_genome_matrix, sparse = TRUE)
  new_jac_genome_matrix <- jaccardMatrix(new_sparse_genome_matrix)
  tmp <- create_new_genome(genome = new_bin_list, sim_matrix = new_jac_genome_matrix)
  found <- tmp$total_solutions
}
new_genome <- tmp$combined_genome

save("file path")

# drawing the progression of the fitness values---- 
# bind multiple dataframes together to plot the change of the fitness functions
my.list2 <- list()
for (j in 1:10) {
  # var <- paste("gen", j, sep="")
  var <- j
  var2 <- get(paste("gen", j, sep = ""))
  my.list2[[var]] <- var2
}
bind_df <- dplyr::bind_rows(my.list2, .id="generations")

# to plot the progress of the objectives over the generations
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
cols <- myColorRamp(c("red", "white"), as.numeric(bind_df$generations))

# plot 3d
fig <- plot_ly(bind_df, x = ~ASR, y = ~SCS, z = ~VET,
               marker = list(color = ~generations, colorscale = cols, showscale = TRUE, size=3))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Average Spearman Rho', autorange="reversed", range = c(0.03, 0.17)),
                                   yaxis = list(title = 'Submatrix Correlation Score', autorange="reversed"),
                                   zaxis = list(title = 'Transposed Virtual Error')),
                      title = "Fitness values over generations",
                      autosize=TRUE,
                      annotations = list(
                        x = 1.13,
                        y = 1.05,
                        text = 'Generations',
                        xref = 'paper',
                        yref = 'paper',
                        showarrow = FALSE
                      ))
fig