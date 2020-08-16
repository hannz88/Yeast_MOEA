# R 3.6.1
# functions to use in nsga2 whitebox

# load dependencies
library(rlist)
library(ecr)

# sort through the parents
get_parents <- function(population, fitness_population){
  # return two list
  # list 1: the population who's fit
  # list 2: the corresponding fitness matix
  pop_idx <- seq_along(population)
  nondom_pop <- doNondominatedSorting(fitness_population)
  max_rank <- ceiling(max(nondom_pop$ranks)/2)
  fitness_bin <- as.numeric(nondom_pop$ranks < max_rank)
  parents_idx <- pop_idx[fitness_bin==1]
  new_population <- population[c(parents_idx)]
  newpop_fitness <- fitness_population[, c(parents_idx)]
  tmp <- list(pop = new_population, fit = newpop_fitness)
  return(tmp)
}

# select parents and crossover
get_gametes <- function(new_population, newpop_fitness, lambda){
  # mates two parents and crossover the genome
  gametes <- list()
  for (i in 1:(lambda/2)) {
    parents <- selSimple(newpop_fitness,2)
    p1_idx <- parents[1]
    p2_idx <- parents[2]
    parents_pair <- c(new_population[p1_idx][1], new_population[p2_idx][1])
    recomb <- recUnifCrossover(parents_pair, p=0.7)
    gametes <- list.append(gametes, recomb[1][[1]])
    gametes <- list.append(gametes, recomb[2][[1]])
  }
  return(gametes)
}

# flip the gametes 
get_offspring <- function(gamete){
  offspring <- list()
  for (i in 1:length(gamete)) {
    offspring <- list.append(offspring, mutBitflip(ind = gamete[[i]], p = 0.1))
  }
  return(offspring)
}
  