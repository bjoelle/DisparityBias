

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

simulation.pipeline <- function(){

  biogeography.stuck.count <- 0
  repeat {
  biogeography.stuck <- FALSE
    #  while #while biogeography loop failed, loop through code. If biogeography succeeds, continue
    ### Step 1: Simulate tree 
    # current assumption: the trait value for each branch (i.e. each species) is the value at the end of the branch - decision made to maximise diffs between species
    # current assumption: bifurcating speciation = each branch is a species
    tr <- TreeSim::sim.bd.taxa(n = tips, 1, birth, death)[[1]]
    #plot(tr)
    
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
    #  biogeography.stuck = FALSE
    while(fossils_in_area1 < L || fossils_in_area1 > H) {

      ## Running the biogeography simulation
      traits.bio <- FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = rate.bio)
      ## Updating the number of fossils
      fossils_in_area1 <- sum(traits.bio == '1')
      iteration.count <- iteration.count + 1
      if (iteration.count >= iteration.limit) {
        biogeography.stuck = TRUE
      }
      if (biogeography.stuck == TRUE) {
        break
      }
    }
    
    if (biogeography.stuck == TRUE){
      biogeography.stuck.count <- biogeography.stuck.count + 1
      if (biogeography.stuck.count >= num_rep*2) {
        stop("Stuck in deep biogeography mud")
      }      
    }
    else
    {
      break
    }

  }
  ### Step 5: Simulate biased sampling
  # associate high and low sampling with biogeographical areas in traits.bio [input]
  # simulate biased sampling on tree [input]
  # output is fossils.bio = fossil taxa and respective ages when sampling is biased
  
  rates <- translate.states(traits.bio, low, high)
  
  fossils.bio.dupl <- FossilSim::sim.fossils.poisson(rates, tree = tr)
  
  # colours
  fcols = sapply((traits.bio[unlist(sapply(fossils.bio.dupl$sp, function(i) which(taxa$sp == i)))]), function(j) if(j == 1) fcol1 else fcol2)
  plot(fossils.bio.dupl, tr, strata = bins, show.strata = TRUE, fossil.col = fcols)
  
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

#### Analysis
#TG: You can then use the function replicate to get some replicates. For example:
my_simulations <- replicate(num_rep, simulation.pipeline(), simplify = FALSE)

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
