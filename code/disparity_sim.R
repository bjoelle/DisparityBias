

setwd("~/Desktop/disparity/DisparityBias/code/")

source("functions_DisaBiss.R")


set.seed(17)


### Setting up variables

# Trees
birth <- 0.1 #birth rate
death <- 0.05 # death rate
tips <- 100 # number of tips in tree
# Traits
trait_num <- 2 # number of traits we are simulating
v <- 0.005 # rate of trait evolution
# Uniform Sampling
rate <- 0.2 # rate of fossilisation
# Biogeography simulation
rate.bio = 0.005 # migration rate 
fossils_in_area1 <- 0 # setting up parameter for checking spatial split
threshold <- 0.45 # threshold for spatial split between areas 0 and 1
iteration.limit <- 100 #number of times loop for generating biogeographic areas can loop
# Biased sampling
low = 0.0015 # sampling rate for fossils in low sampling area
high = 0.5 # sampling rate for fossils in high sampling area
# Time binning
bins <- 3 # number of time bins

simulate_one_set <- function(){
  
  ### Step 1: Simulate tree 
  # current assumption: the trait value for each branch (i.e. each species) is the value at the end of the branch - decision made to maximise diffs between species
  # current assumption: bifurcating speciation = each branch is a species
  tr <- TreeSim::sim.bd.taxa(n = tips, 1, birth, death)[[1]]
  plot(tr)
  
  taxa <- FossilSim::sim.taxonomy(tr, beta = 1) # how you define morphotaxa with respect to the tree
  
  
  ### Step 2: Simulate "true" disparity
  traits <- generate.traits(taxa, trait_num, tr, v)
  
  ### Step 3: Simulate constant rate of preservation
  
  fossils.uni.dupl <- FossilSim::sim.fossils.poisson(rate = rate, tree = tr)
  plot(fossils.uni.dupl, tr, strata = bins, show.strata = TRUE)
 
  ### Step 4: Simulate biogeography on tree
  # assumption: approach assumes migration does not influence tree shape
  # calculate threshold values for number of taxa in each geographic area
  # loop to simulate biogeography as a binary character under the Mk model
  # loop keeps track of number of attempts, exits if iteration.limit is reached
  # inputs tree and migration rate
  # outputs traits.bio, object of type double
  
  #TO DO: integrate resetting iteration.count before 'if'
  
  number_of_tips <- length(traits$sp)
  L <- round(sum(threshold*number_of_tips))
  H <- sum(number_of_tips-L)
  iteration.count <- 0 #always set at 0, resetting it
  while(fossils_in_area1 < L || fossils_in_area1 > H) {
    if (iteration.count >= iteration.limit) {
      stop("Failed to converge on a suitable geographical distribution")
    }
    ## Running the biogeography simulation
    traits.bio <- FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = rate.bio)
    ## Updating the number of fossils
    fossils_in_area1 <- sum(traits.bio == '1')
    iteration.count <- iteration.count + 1
  }
  
  ### Step 5: Simulate biased sampling
  # associate high and low sampling with biogeographical areas in traits.bio [input]
  # simulate biased sampling on tree [input]
  # output is fossils.bio = fossil taxa and respective ages when sampling is biased
  
  rates <- translate.states(traits.bio, low, high)
  
  fossils.bio.dupl <- FossilSim::sim.fossils.poisson(rates, tree = tr)
  plot(fossils.bio.dupl, tr, strata = bins, show.strata = TRUE)
  
  ### Step 6: Bin fossils and match traits with species & bins
  # assumption: no extant samples simulated or sampled, although some fossil species may be extant 
  # calculate bin max/min ages based on tree [input] and number of bins
  # 
  
  max.age <- FossilSim::tree.max(tr)
  int.ages <- seq(0, max.age, length = bins + 1)
  
  ###### run Joelle's function
  boop <- bin.taxa(taxa, 3, max.age)
  all.binned <- FossilSim::sim.interval.ages(boop, max.age = max.age, strata = bins, use.species.ages = FALSE)
  
  # bin fossils for unbiased sampling set
  fossils.binned <- FossilSim::sim.interval.ages(fossils.uni.dupl, tr, max.age = max.age, strata = bins, use.species.ages = FALSE)
  # bin fossils for biased sampling set
  fossils.bio.binned <- FossilSim::sim.interval.ages(fossils.bio.dupl, tr, max.age = max.age, strata = bins, use.species.ages = FALSE)
  
  
  
  bias <- int.assign(fossils.bio.binned, int.ages) ##
  uni <- int.assign(fossils.binned, int.ages)
  all <- int.assign(all.binned, int.ages)
  
  #return(list(bias = bias, uni = uni, all = all,  
  #            traits = traits, traits.bio = traits.bio))
  
  ## Rename some variables
  my_geography <- traits.bio
  #my_geography <- data.frame(traits.bio, taxa$sp) ## change to adding bins not taxa numbers
  my_trait_space <- traits[, c("trait1", "trait2")]
  
  uni_sample <- subset(uni$sp, uni$int == "2")
  uni_sample <- uni_sample[!duplicated(uni_sample)] #removing duplicates
  
  bias_sample <- subset(bias$sp, bias$int == "2")
  bias_sample <- bias_sample[!duplicated(bias_sample)]
  
  ## Creating the group vector for dispRity
  my_groups <- list(## All the species
    "all_species" = subset(all$sp, all$int == "2"),
    ## All species in location 1
    "area_0" = subset(all$sp, all$int == "2")[(subset(all$sp, all$int == "2") %in% which(my_geography == 0))],
    ## All species in location 2
    "area_1" = subset(all$sp, all$int == "2")[(subset(all$sp, all$int == "2") %in% which(my_geography == 1))],
    ## The uniform sampled group
    "uni_sample" = uni_sample,
    ## The biased sampled group
    "bias_sample" = bias_sample,
    ## unif species sampling
    "uni_species" = sample(subset(all$sp, all$int == "2"), 20)
    )
  
  ## Creating a dispRity object that contains the trait space and the groups
  my_groupings <- custom.subsets(data = as.matrix(my_trait_space),
                                 group = my_groups)
  #TG: ignore the warning (or read it to know what it just did ;) - but nothing bad happening here)
  
  return(my_groupings)
  
}

tmp = simulate_one_set()

# create a new data.frame for disparity analyses based on sampled species in each bin, puts data in format for disparity analysis
#disparity.df <- function(traits, fossils, interval.ages){
#  if(identical(fossils$hmin, fossils$hmax))
#    stop("fossils must be binned!")
#  disp <- data.frame()
#  for(i in 1:(bins - 1)){
#    hmin <- interval.ages[i]
#    hmax <- interval.ages[i+1]
#    tmp <- subset(fossils, hmin == interval.ages[i])
#    tmp <- unique(tmp)
#    for(j in tmp$sp){
#      tmp2 <- data.frame(bin = i, sp = j, bin.midpoint = mean(c(hmin, hmax)))
#      for(k in 1:trait_num){
#        t <- traits[which(traits$sp == j),][paste0("trait",k)][[1]]
#        tmp2 <- cbind(tmp2, data.frame(tc = t))
#        colnames(tmp2)[ncol(tmp2)] <-  paste0("trait",k)
#      }
#      disp <- rbind(disp, tmp2) 
#    }
#  }
#  disp
#}

# disp <- disparity.df(traits, fossils.binned, int.ages) ## ? sampled fossils with uniform sampling
# disp.bio <- disparity.df(traits, fossils.bio.binned, int.ages) ## ? sampled fossils with biased sampling
# 
# h1 <- hist(traits$trait1)
# h2 <- hist(disp$trait1, breaks = h1$breaks)
# h3 <- hist(disp.bio$trait1, breaks = h1$breaks)
# 
# plot( h1, col=rgb(0,0,1,1/4), xlab = "Trait values", main = NULL)  # first histogram
# plot( h2, col=rgb(1,0,0,1/4), add=T)  # second
# plot( h3, col=rgb(0,0.8,0.6,1/4), add=T)  # third
# 
# legend(x = "topright", legend = c("True disparity", "Uniform sampling","Biased sampling"),
#        fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), rgb(0,.5,.5,1/4)) )
# 
# plot( (((traits$start - traits$end) / 2) + traits$start), traits$trait1, col = rgb(0,0,1,1/4), pch = 19)
# points(disp$bin.midpoint, disp$trait1, col = rgb(1,0,0,1/4), pch = 19)
# points(disp.bio$bin.midpoint, disp.bio$trait1, col = rgb(0,.5,.5,1/4), pch = 19)


### Step 7: Compute disparity metrics

# # compute metric - centroids
# library(dispRity)
# disparity_centr <- dispRity::dispRity.per.group(ordinated_all,
#                                      list(trueO = c(1:true.mx), uni = c(uni.mn:uni.mx), bias = c(bias.mn:bias.mx)),
#   metric = centroids) #variance
# # metric = variances)
# disparity_centr
# 
# plot(disparity_centr)
# 
# #compute metric - sum of variances
# ####c(mean,variance))
# disparity_var <- dispRity::dispRity.per.group(ordinated_all,
#                                                list(trueO = c(1:true.mx), uni = c(uni.mn:uni.mx), bias = c(bias.mn:bias.mx)),
#                                                metric = variances) #variance
# # metric = variances)
# disparity_var
# plot(disparity_var)
# plot(disparity_centr, type = "preview")


my_groupings = tmp



## Calculating disparity on these groups
disparity_sum_var <- dispRity(my_groupings, metric = c(sum, variances))

## Hop
plot(disparity_sum_var)
#TG: note that these are now point estimates (one disparity value per group) hence the absence of variance. You'll get the variance from replicating the simulations (see pseudo code below).
#TG: to just get the values you can use the function get.disparity (see example in the pseudo code below)
#get.disparity(disparity_sum_var)
extract.dispRity(disparity_sum_var)

scatter_data <- cbind(my_trait_space, my_geography)
ggplot(my_trait_space, aes(x=trait1, y=trait2, color=my_geography)) + geom_point()


## Write pseudo function for the simulation pipeline
simulation.pipeline <- function(birth, death, tips, <other_magic_numbers>) {
  ## Simulating the data
  <the pipeline that simulates the species, the traits and the sampling bias>
  ## Preparing the data
  <the code to create the trait space to be fed to dispRity, like the code from above with my_groups and custom.subsets>
  ## The output
  return(<the bit to be returned, again, I suggest outputing the my_groupings variable from the example above>)
}

#TG: You can then use the function replicate to get some replicates. For example:
my_5_simulations <- replicate(5, simulation.pipeline(birth = X, death = Y, tips = Z, <other_magic_numbers>))

#TG: You can then measure the disparity on the output using lapply (applying a function to a list)
my_sum_variances <- lapply(my_5_simulations, dispRity, metric = c(sum, variances))

#TG: You can then extract the disparity values (the point estimates explained above)
my_point_estimates <- lapply(my_sum_variances, get.disparity)

##TG: And you can combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
my_results_table <- do.call(rbind, my_point_estimates)
boxplot(my_results_table)
