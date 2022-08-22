

setwd("~/Desktop/disparity/DisparityBias/code/")

source("functions_DisaBiss.R")


set.seed(17)


### Setting up variables

# Trees
birth <- 0.1 # birth rate
death <- 0.05 # death rate
tips <- 10 # number of tips in tree
# Traits
trait_num <- 2 # number of traits we are simulating
v <- 0.005 # rate of trait evolution
# Uniform Sampling
rate <- 0.2 # rate of fossilisation
# Biogeography simulation
rate.bio = 0.05 # migration rate 
fossils_in_area1 <- 0 # setting up parameter for checking spatial split
threshold <- 0.45 # threshold for spatial split between areas 0 and 1
iteration.limit <- 100 #number of times loop for generating biogeographic areas can loop
# Biased sampling
low = 0.0015 # sampling rate for fossils in low sampling area
high = 0.5 # sampling rate for fossils in high sampling area
# Time binning
bins <- 3 # number of time bins
#Colours for fossils in tree plots
fcol1 <- "#FF6EB4"
fcol2 <- "#C0FF3E"

num_rep <- 5


#### Analysis
#TG: You can then use the function replicate to get some replicates. For example:
my_simulations <- replicate(num_rep, simulation.pipeline(birth, death, tips, trait_num, v, rate, rate.bio, fossils_in_area1, threshold, iteration.limit, low, high, bins, fcol1, fcol2), simplify = FALSE)

#TG: You can then measure the disparity on the output using lapply (applying a function to a list)
my_sum_variances <- lapply(my_simulations, dispRity, metric = c(sum, variances))

#TG: You can then extract the disparity values (the point estimates explained above)
my_point_estimates <- lapply(my_sum_variances, extract.dispRity)

##TG: And you can combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
my_results_table <- do.call(rbind, my_point_estimates)

columns <- c("values", "sampling")
result <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(result) <- columns

for (i in 1:ncol(my_results_table)){
  for (j in my_results_table[,i]){
    current.row = c(j,colnames(my_results_table)[i])
    result[nrow(result) + 1,] <- current.row
  }
}

result$sampling <- as.factor(result$sampling)
result$values <- as.numeric(result$values)
p <- ggplot(result, aes(x=sampling, y=values)) +
  geom_boxplot()
p
