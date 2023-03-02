
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
migration.rate_multi <- c(0.001, 0.003, 0.01) # migration rate*
threshold <- 0.45 # threshold for spatial split between areas 0 and 1*

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

num.rep <- 10

sims = TRUE
analysis = TRUE

# a place to store output
if(!dir.exists(outdir)) dir.create(outdir)

### Simulations

# checks which variable contains multiple values and loops through them
var = strsplit( ls(pat = "multi"), "_" )[[1]][1]
vals = eval( parse( text = ls(pat = "multi") ) )

for(i in vals){
  
  # assign 
  assign( var, i ) 
  if(sims){
    #TODO: add morphospace plots to the output
    simulations <- lapply(1:num.rep, function(x){simulation.pipeline(birth, death, tips, trait.num, trait.evol.rate, fossilisation.rate, migration.rate, fossils.in.area1, threshold, iteration.limit, low.sampling, high.sampling, bins, fossil.colour1, fossil.colour2, x, var, i)})
    save(simulations, file = paste0(outdir, "data_", var, "_", i, "_", ".RData")) #TODO: need a naming convention for different simulation conditions
  } else {
    load(file = paste0(outdir, "data_", var, "_", i, "_", ".RData"))
  }
  
  #TODO: print this to file
  # Check if enough samples present in subsamples
  for (j in 1:num.rep){
    if(lengths(simulations[[j]]$subsets$area_0) < 40 || lengths(simulations[[j]]$subsets$area_1) < 40) {
      print(paste("Too few fossils in run", j))
    }
    else print(paste("All good", j))
  }
  
  if(analysis){
    ### Disparity Analysis - these functions return plots
    sumv <- disparity.analysis(simulations, analysis = "sum of variances")
    mpd <- disparity.analysis(simulations, analysis = "pairwise distance")
    mcd <- disparity.analysis(simulations, analysis = "centroids")
    
    assign(paste0("sumv_", var, "_", i), sumv)
    assign(paste0("mpd_", var, "_", i), mpd)
    assign(paste0("mcd_", var, "_", i), mcd)
  }
}

wd = 10
ht = 3

pdf(file = paste0("sumv_", var, "_.pdf"), width = wd, height = ht)
par(mfcol=c(1, 3))
sumv_migration.rate_0.001
sumv_migration.rate_0.003
sumv_migration.rate_0.01
dev.off()

pdf(file = paste0("mpd_", var, "_.pdf"), width = wd, height = ht)
par(mfcol=c(1, 3))
mpd_migration.rate_0.001
mpd_migration.rate_0.003
mpd_migration.rate_0.01
dev.off()

pdf(file = paste0("mcd_", var, "_.pdf"), width = wd, height = ht)
par(mfcol=c(1, 3))
mcd_migration.rate_0.001
mcd_migration.rate_0.003
mcd_migration.rate_0.01
dev.off()



