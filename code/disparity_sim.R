

setwd("~/Desktop/DisparityBias/code/")

source("functions_DisaBiss.R")


set.seed(20)


### Setting up variables

# Trees
birth <- 0.1 # birth rate
death <- 0.05 # death rate
tips <- 100 # number of tips in tree
# Traits
trait.num <- 2 # number of traits we are simulating
trait.evol.rate <- 0.005 # rate of trait evolution
# Uniform Sampling
fossilisation.rate <- 0.1 # rate of fossilisation
# Biogeography simulation
migration.rate <- 0.05 # migration rate 
fossils.in.area1 <- 0 # setting up parameter for checking spatial split
threshold <- 0.45 # threshold for spatial split between areas 0 and 1
iteration.limit <- 100 #number of times loop for generating biogeographic areas can loop
# Biased sampling
low.sampling <- 0.0005 # sampling rate for fossils in low sampling area
high.sampling <- 0.1 # sampling rate for fossils in high sampling area
# Time binning
bins <- 3 # number of time bins
#Colours for fossils in tree plots
fossil.colour1 <- "#FF6EB4"
fossil.colour2 <- "#C0FF3E"

num.rep <- 2

#Running the simulations
simulations <- replicate(num.rep, simulation.pipeline(birth, death, tips, trait.num, trait.evol.rate, fossilisation.rate, migration.rate, fossils.in.area1, threshold, iteration.limit, low.sampling, high.sampling, bins, fossil.colour1, fossil.colour2), simplify = FALSE)



#### Analysis
#Measure the disparity on the output using lapply (applying a function to a list)
sum.variances <- lapply(simulations, dispRity, metric = c(sum, variances))

#Extract the disparity values (the point estimates explained above)
point.estimates <- lapply(sum.variances, extract.dispRity)

##Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
results.table <- do.call(rbind, point.estimates)

columns <- c("values", "sampling")
result <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(result) <- columns

for (i in 1:ncol(results.table)){
  for (j in results.table[,i]){
    current.row = c(j,colnames(results.table)[i])
    result[nrow(result) + 1,] <- current.row
  }
}

result$sampling <- as.factor(result$sampling)
result$values <- as.numeric(result$values)
p <- ggplot(result, aes(x=sampling, y=values)) +
  geom_boxplot()
p


