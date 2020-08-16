# R 3.6.1
# nsga2 in exterior using newly improved fitness functions with cpp
# 15th July, 2020

# load dependencies----
if(!require(rlist)){
  install.packages("rlist")
  library(rlist)
}
if(!require(ecr)){
  install.packages("ecr")
  library(ecr)
}
if(!require(Rcpp)){
  install.packages("Rcpp")
  library(Rcpp)
}

# make sure you run R CMD INSTALL yeastfitness in the local terminal where the 
# yeastfitness file is stored first 
require(yeastfitness)

# load the file that has the functions to use nsga2 as white box----
source("file path for NSGA2_subfunctions.R")

# load files----
data <- read.csv("file path for dataset", 
                     sep=",", header = TRUE, row.names = 1)
data <- na.omit(data)
data <- data[-5]  # to remove cdc28

# combine all the functions into a fitness functions----
fitness <- function(chromosome){
  row_vect <- chromosome[1:nrow(data)]
  col_vect <- chromosome[-(1:nrow(data))]
  if ((sum(row_vect) <=1) | (sum(col_vect) <=1)) {
    result <- c(0, 100, 100)
  } else {
    mx <- data[c(row_vect==1), c(col_vect==1)]  # bicluster
    matrix_mx <- as.matrix(mx)
    f1 <- yeastfitness::get_ASR_cpp(matrix_mx)
    f2 <- yeastfitness::get_SCS_cpp(matrix_mx)
    f3 <- yeastfitness::get_VET_cpp(matrix_mx)
    result <- c(f1, f2, f3)
  }
  return(result)
}


# running nsga2 as a white box----
mu = 200
lambda = 200
max_iter <- 15000
var_num <- nrow(data) + ncol(data)

# init control
control <- initECRControl(fitness.fun = fitness, n.objectives = 3, minimize = c(FALSE, TRUE, TRUE))
control = registerECROperator(control, "selectForSurvival", selNondom)  # use this for Non-dominated sorting

# initiate population
init_population <- genBin(mu, n.dim = var_num)  # n.dim = the length of the chromosome
fitness_pop <- evaluateFitness(control, init_population)


is.error <- function(x) inherits(x, "try-error")  # try catch
# evolutionary method   # this method will print out an RData file for every generation
# useful if you need to plot the progression
start_time <- Sys.time()  # start time
for (i in seq_len(max_iter)) {  # only works if you're testing 1 gen
  retries <- 0
  tmp_name <- paste0("gen",i)
  fitness_pop <- evaluateFitness(control, init_population)
  parents <- get_parents(population=init_population, fitness_population = fitness_pop)  # select parents by selecting the fittest half
  parents_genome <- parents$pop 
  parents_fitness <- parents$fit
  gametes_genome <- get_gametes(new_population = parents_genome, newpop_fitness = parents_fitness)  # get gametes from the parents
  while (retries < 3) {  # some times the algorithm gets stuck 
    offspring <- get_offspring(gametes_genome) # get offspring 
    fitness_offspring <- evaluateFitness(control, offspring)
    pool <- try(replaceMuPlusLambda(control, init_population, offspring, fitness_pop, fitness_offspring))
    FAILED <- is.error(pool)
      if (FAILED == TRUE) {
      retries = retries + 1
      rdata_name <- paste0("file path/output_", tmp_name, "_error.RData")
      save.image(rdata_name)
    } else {
      break
    }
  }
  init_population <- pool$population
  fitness_pop <- pool$fitness
  tmp_df <- as.data.frame(t(fitness_pop))  # save the fitness values of solution from each generation
  colnames(tmp_df) <- c("ASR", "SCS", "VET")  # assign column names to the dataframe
  assign(tmp_name, tmp_df)  # save the dataframe to a different name
  end_time <- Sys.time()  # end time
  rdata_name <- paste0("file path/output_", tmp_name, ".RData")
  save.image(rdata_name)
} 
