
### functions and libraries
source("functions_DisaBiss.R")

library(dispRity)
library(FossilSim)
library(TreeSim)
library(ggplot2)

set.seed(23)

outdir="/Users/warnock/Documents/files/work/research/projects/DisparityBias/output/"

### Setting up variables

# Trees
birth <- 0.1 # birth rate
death <- 0.075 # death rate
tips <- 200 # number of tips in tree

# Traits
trait.num <- 2 # number of traits we are simulating
trait.evol.rate <- 0.01 # rate of trait evolution* 0.03 - fine

# Uniform Sampling
fossilisation.rate <- 0.05 # rate of fossilisation

# Biogeography simulation
migration.rate <- 0.003 # migration rate*
threshold <- 0.45 # threshold for spatial split between areas 0 and 1

fossils.in.area1 <- 0 # setting up parameter for checking spatial split
iteration.limit <- 100 #number of times loop for generating biogeographic areas can loop

# Biased sampling
low.sampling <- 0.01 # sampling rate for fossils in low sampling area*
high.sampling <- 0.1 # sampling rate for fossils in high sampling area

# Time binning
bins <- 3 # number of time bins

#Colours for fossils in tree plots
fossil.colour1 <- "#5AA8C5"
fossil.colour2 <- "#F8D754"

num.rep <- 20

### Running the simulations RW: switched to lapply so I can use the iteration number within the function
simulations <- lapply(1:num.rep, function(x){simulation.pipeline(birth, death, tips, trait.num, trait.evol.rate, fossilisation.rate, migration.rate, fossils.in.area1, threshold, iteration.limit, low.sampling, high.sampling, bins, fossil.colour1, fossil.colour2, x)})
if(!dir.exists(outdir)) dir.create(outdir)
save(simulations, file = paste0(outdir, "data.RData")) #TODO: need a naming convention for different simulation conditions

# Check if enough samples present in subsamples
for (i in 1:num.rep){
  if(lengths(simulations[[i]]$subsets$area_0) < 40 || lengths(simulations[[i]]$subsets$area_1) < 40) {
    print(paste("Too few fossils in run", i))
  }
  else print(paste("All good", i))
}

### Analysis
sumv <- disparity.analysis(simulations, analysis = "sum of variances")
mpd <- disparity.analysis(simulations, analysis = "pairwise distance")
mcd <- disparity.analysis(simulations, analysis = "centroids")

sumv
mpd
mcd

# sum.variances <- lapply(simulations, dispRity, metric = c(sum, variances))
# 
# ####
# # Extract the disparity values (the point estimates explained above)
# point.estimates.sumv <- lapply(sum.variances, get.disparity)
# 
# ## Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
# results.table.sumv <- do.call(rbind, point.estimates.sumv)
# 
# columns <- c("values", "sampling")
# result.sumv <- data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(result.sumv) <- columns
# 
# for (i in 1:ncol(results.table.sumv)){
#   for (j in results.table.sumv[,i]){
#     current.row = c(j,colnames(results.table.sumv)[i])
#     result.sumv[nrow(result.sumv) + 1,] <- current.row
#   }
# }
# 
# 
# # observed_disparity - mean(null_disparity)
# null_disparity <- mean(unlist(results.table.sumv[,1]))
# 
# # sampling regime
# result.sumv$sampling <- as.factor(result.sumv$sampling)
# # disparity
# result.sumv$values <- as.numeric(result.sumv$values) - null_disparity
# 
# # plots
# p <- ggplot(result.sumv, aes(x=sampling, y=values)) +
#   labs(title = "Sum of Variances") +
#   geom_boxplot()
# p
# 
# 
# #### Analysis -- Median pairwise distances
# #Measure the disparity on the output using lapply (applying a function to a list)
# med.pair.distances <- lapply(simulations, dispRity, metric = c(median, pairwise.dist))
# 
# #Extract the disparity values (the point estimates explained above)
# point.estimates.mpd <- lapply(med.pair.distances, get.disparity)
# 
# ##Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
# results.table.mpd <- do.call(rbind, point.estimates.mpd)
# 
# columns <- c("values", "sampling")
# result.mpd <- data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(result.mpd) <- columns
# 
# for (i in 1:ncol(results.table.mpd)){
#   for (j in results.table.mpd[,i]){
#     current.row = c(j,colnames(results.table.mpd)[i])
#     result.mpd[nrow(result.mpd) + 1,] <- current.row
#   }
# }
# 
# result.mpd$sampling <- as.factor(result.mpd$sampling)
# result.mpd$values <- as.numeric(result.mpd$values)
# q <- ggplot(result.mpd, aes(x=sampling, y=values)) +
#   labs(title = "Median pairwise distances") +
#   geom_boxplot()
# q
# 
# 
# 
# #### Analysis -- Median distance from centroid (centroid = origin)
# #Measure the disparity on the output using lapply (applying a function to a list)
# med.centroid.distances <- lapply(simulations, dispRity, metric = c(median, centroids), centroid = 0)
# 
# #Extract the disparity values (the point estimates explained above)
# point.estimates.mcd <- lapply(med.centroid.distances, get.disparity)
# 
# ##Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
# results.table.mcd <- do.call(rbind, point.estimates.mcd)
# 
# columns <- c("values", "sampling")
# result.mcd <- data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(result.mcd) <- columns
# 
# for (i in 1:ncol(results.table.mcd)){
#   for (j in results.table.mcd[,i]){
#     current.row = c(j,colnames(results.table.mcd)[i])
#     result.mcd[nrow(result.mcd) + 1,] <- current.row
#   }
# }
# 
# result.mcd$sampling <- as.factor(result.mcd$sampling)
# result.mcd$values <- as.numeric(result.mcd$values)
# c <- ggplot(result.mcd, aes(x=sampling, y=values)) +
#   labs(title = "Median distance from centroids") +
#   geom_boxplot()
# c
