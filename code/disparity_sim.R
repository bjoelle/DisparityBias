

setwd("~/Desktop/disparity/DisparityBias/code/")

source("functions_DisaBiss.R")


set.seed(20)


### Setting up variables

# Trees
birth <- 0.1 # birth rate
death <- 0.05 # death rate
tips <- 1000 # number of tips in tree
# Traits
trait.num <- 2 # number of traits we are simulating
trait.evol.rate <- 0.001 # rate of trait evolution
# Uniform Sampling
fossilisation.rate <- 0.05 # rate of fossilisation
# Biogeography simulation
migration.rate <- 0.03 # migration rate 
fossils.in.area1 <- 0 # setting up parameter for checking spatial split
threshold <- 0.45 # threshold for spatial split between areas 0 and 1
iteration.limit <- 100 #number of times loop for generating biogeographic areas can loop
# Biased sampling
low.sampling <- 0.0003 # sampling rate for fossils in low sampling area
high.sampling <- 0.08 # sampling rate for fossils in high sampling area
# Time binning
bins <- 3 # number of time bins
#Colours for fossils in tree plots
fossil.colour1 <- "#FF6EB4"
fossil.colour2 <- "#C0FF3E"

num.rep <- 5

### Running the simulations
simulations <- replicate(num.rep, simulation.pipeline(birth, death, tips, trait.num, trait.evol.rate, fossilisation.rate, migration.rate, fossils.in.area1, threshold, iteration.limit, low.sampling, high.sampling, bins, fossil.colour1, fossil.colour2), simplify = FALSE)

# Check if enough samples present in subsamples
for (i in 1:num.rep){
  if (lengths(simulations[[i]]$subsets$area_0) < 40 || lengths(simulations[[i]]$subsets$area_1) < 40) {
    print("Too few fossils in run")
  }
  else print("All good")
}




### Scatter plot

all_species <- simulations[[1]]$subsets$all_species
area_0 <- simulations[[1]]$subsets$area_0
area_1 <- simulations[[1]]$subsets$area_1
bias_0 <- simulations[[1]]$subsets$bias_0_sample
bias_1 <- simulations[[1]]$subsets$bias_1_sample

traits <- simulations[[1]]$matrix

#plot all traits
traits.df <- as.data.frame(traits)
a <-ggplot(traits.df, aes(x=trait1, y=trait2)) + geom_point()

#add column to data frame
len <- nrow(traits.df)
geog <- rep(5,len)
traits.df <- cbind(traits.df, geog = geog) #prepopulate with 0s

#loop through dataframe, if row should be area 1 ten change column value

for (j in 1:nrow(traits.df)){
  traits.df$geog <- (traits.df$geog + 1)
}



bob <- sapply(geog, function(j) if )

%in% which(rownames(traits.df) == bias_0)




#### Analysis -- Sum of variances
# Measure the disparity on the output using lapply (applying a function to a list)
sum.variances <- lapply(simulations, dispRity, metric = c(sum, variances))

# Extract the disparity values (the point estimates explained above)
point.estimates.sumv <- lapply(sum.variances, extract.dispRity)

## Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
results.table.sumv <- do.call(rbind, point.estimates.sumv)

columns <- c("values", "sampling")
result.sumv <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(result.sumv) <- columns

for (i in 1:ncol(results.table.sumv)){
  for (j in results.table.sumv[,i]){
    current.row = c(j,colnames(results.table.sumv)[i])
    result.sumv[nrow(result.sumv) + 1,] <- current.row
  }
}

result.sumv$sampling <- as.factor(result.sumv$sampling)
result.sumv$values <- as.numeric(result.sumv$values)
p <- ggplot(result.sumv, aes(x=sampling, y=values)) +
  labs(title = "Sum of Variances") +
  geom_boxplot()
p




#### Analysis -- Median pairwise distances
#Measure the disparity on the output using lapply (applying a function to a list)
med.pair.distances <- lapply(simulations, dispRity, metric = c(median, pairwise.dist))

#Extract the disparity values (the point estimates explained above)
point.estimates.mpd <- lapply(med.pair.distances, extract.dispRity)

##Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
results.table.mpd <- do.call(rbind, point.estimates.mpd)

columns <- c("values", "sampling")
result.mpd <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(result.mpd) <- columns

for (i in 1:ncol(results.table.mpd)){
  for (j in results.table.mpd[,i]){
    current.row = c(j,colnames(results.table.mpd)[i])
    result.mpd[nrow(result.mpd) + 1,] <- current.row
  }
}

result.mpd$sampling <- as.factor(result.mpd$sampling)
result.mpd$values <- as.numeric(result.mpd$values)
q <- ggplot(result.mpd, aes(x=sampling, y=values)) +
  labs(title = "Median pairwise distances") +
  geom_boxplot()
q



#### Analysis -- Median distance from centroid (centroid = origin)
#Measure the disparity on the output using lapply (applying a function to a list)
med.centroid.distances <- lapply(simulations, dispRity, metric = c(median, centroids), centroid = 0)

#Extract the disparity values (the point estimates explained above)
point.estimates.mcd <- lapply(med.centroid.distances, extract.dispRity)

##Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
results.table.mcd <- do.call(rbind, point.estimates.mcd)

columns <- c("values", "sampling")
result.mcd <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(result.mcd) <- columns

for (i in 1:ncol(results.table.mcd)){
  for (j in results.table.mcd[,i]){
    current.row = c(j,colnames(results.table.mcd)[i])
    result.mcd[nrow(result.mcd) + 1,] <- current.row
  }
}

result.mcd$sampling <- as.factor(result.mcd$sampling)
result.mcd$values <- as.numeric(result.mcd$values)
c <- ggplot(result.mcd, aes(x=sampling, y=values)) +
  labs(title = "Median distance from centroids") +
  geom_boxplot()
c
