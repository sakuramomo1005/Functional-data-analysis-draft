########################################
# hand-write genetic algorithm 
# date: 2019-3-4
########################################



################ Functions ###################
# just run the follwing functions

# generation an initial generation
init_population = function(n_population, DNA_size){
  population = matrix(round(runif(n_population * DNA_size,0,1)), n_population, DNA_size) 
  return(population)
}
#population =init_population(10,8)

# change DNA from binary decimal to decimal system.
transformDNA = function(population){
  return((population %*% as.matrix(2^(DNA_size:1 - 1),DNA_size,1) / 2^(DNA_size) - 0.5) * 
           (x_bounder[2] - x_bounder[1]) + 0.5 * (x_bounder[2] + x_bounder[1]) )
}
#transformDNA(population)

# calculate the fitness for each individual 
fitness = function(population){
  transform_population = transformDNA(population)
  fitness_score = c()
  for(ii in 1:length(transform_population)){
    #print(ii)
    fitness_score = c(fitness_score, f(transform_population[ii]))
  }
  fitness_score = matrix(fitness_score, length(fitness_score),1)
  return(fitness_score - min(fitness_score)) 
  # since we need to sample based on the fitness, make it positive
}
fitness_score = fitness(population)

# sample from the population, individuals with high fitness are more likely to be chosen
selection = function(population, fitness_score){
  fitness_score = fitness_score + 1e-4  # ensure that fitness score larger than 0
  idx = sample(1:n_population, n_population, replace = TRUE, p = fitness_score/sum(fitness_score))
  return(population[idx,])
}
selection(population, fitness_score)

# mate, crossover
create_child = function(parent, population){
  if(runif(1,0,1) < cross_rate){
    index = sample(1:n_population,1)
    cross_points = which(sample(c(0,1), DNA_size, replace = TRUE) == 1)
    parent[cross_points] = population[index,cross_points]
  }
  return(parent)
}
parent = population[1,]
create_child(parent, population)

# gene mutation
mutate_child = function(child){
  for(i in 1:DNA_size){
    if(runif(1,0,1) < mutate_rate){
      child[i] = 1 - child[i]
    }else{
      child[i] = child[i]
    }
  }
  return(child)
}
mutate_child(parent)

# evolution
evolution = function(n_population, DNA_size, x_bounder, cross_rate, mutate_rate, n_iterations = 1000){
  population = init_population(n_population, DNA_size)
  for(i in 1:n_iterations){
    #fitness_score = fitness(population)
    # for(jj in 1:n_population){
    #   fitness_score = c(fitness_score, fitness(population[jj,]))
    # }
    fitness_score = fitness(population)
    best_person = population[which(fitness_score == max(fitness_score))[1],]
    
    if(i %% 100 == 1){
      print(paste(i,'th evolution best score is ', f(transformDNA(best_person))))
    }
    
    population = selection(population, fitness_score)
    population_copy = population
    new_population = population
    
    for(n_parent in 1:n_population){
      #print(n_parent)
      parent = population[n_parent,]
      child = create_child(parent, population_copy)
      child = mutate_child(child)
      new_population[n_parent,] = child
    }
    population = new_population
  }
  
  
  #print(best_person)
  print(paste('max point is',transformDNA(best_person)))
  x = transformDNA(best_person) 
  y = f(x)
  return(list(x = x, y = y))
}


################ Examples ###################


n_population = 300
DNA_size = 8
x_bounder = c(-3,4)
cross_rate = 0.9 # the common value in this algorithm
mutate_rate = 0.1 # the common value in this algorithm

f = function(x){return(sin(x) + cos(x))}
result = evolution(n_population, DNA_size, x_bounder, cross_rate, mutate_rate, n_iterations = 200)
test_value = f(seq(-3,4,0.001))
plot(seq(-3,4,0.001), test_value, type = 'l')
points(result$x, result$y, col = 'red', pch = 18)

# we can see the max value is found.