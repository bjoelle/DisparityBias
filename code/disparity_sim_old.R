
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
