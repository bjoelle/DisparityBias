
# generate new file for storing traits with taxa in it already [input]
# simulate trait.num number of traits and append to traits file [output]

generate.traits <- function(taxa, trait.num, tr, trait.evol.rate){
  traits <- taxa
  for(i in 1:trait.num){
    tmp <- FossilSim::sim.trait.values(init = 5, tree = tr, model = "BM", v = trait.evol.rate, min.value = 0)
    traits <- cbind(traits, tmp)
    colnames(traits)[ncol(traits)] <- paste0("trait",i)
  }
  return(traits)
}




# associate high and low sampling with biogeographical areas in fossil.biogeographic.area [input]
translate.states <- function(fossil.biogeographic.area, low.sampling, high.sampling) sapply(fossil.biogeographic.area, function(t) if(t == 1) low.sampling else high.sampling)


#turns all taxa into a format useable by FossilSim so that time binning can occur
bin.taxa = function(taxa, nbins, max.age) {
  if(nbins%%1 != 0 || nbins == 0 || nbins < 0) {
    stop("Number of bins must be a positive integer, check nbins")
  }
  
  bin.ages = seq(0, max.age, max.age/nbins)
  bin.min = bin.ages[-length(bin.ages)]
  bin.max = bin.ages[-1]
  
  fs = data.frame()
  for(tax in 1:nrow(taxa)) {
    overlap = which(bin.min <= taxa$start[tax] & bin.max >= taxa$end[tax])
    if(length(overlap) == 0) {
      warning("Taxa overlaps with no bins, check max.age")
      next
    }
    
    for(bin.index in overlap) {
      min.age = max(bin.min[bin.index], taxa$end[tax])
      max.age = min(bin.max[bin.index], taxa$start[tax])
      mid.age = (min.age + max.age) / 2
      fs = rbind(fs, data.frame(sp = taxa$sp[tax], edge = taxa$edge[tax], hmin = mid.age, hmax = mid.age))
    }
  }
  
  FossilSim::fossils(fs)
}


### function to turn sim.interval.ages into defined/numbered time bins
int.assign <- function(fossils, ints){
  if(identical(fossils$hmin, fossils$hmax))
    stop("fossils must be binned!")
  fossils$int <- NA
  for(i in 1:(length(ints) - 1)){
    if(any(fossils$hmin == ints[i])){
      fossils[which(fossils$hmin == ints[i]),]$int = i
    }
  }
  fossils
}



simulation.pipeline <- function(birth, death, tips, trait.num, trait.evol.rate, fossilisation.rate, migration.rate, fossils.in.area1, threshold, iteration.limit, low.sampling, high.sampling, bins, fossil.colour1, fossil.colour2){
  
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
    traits <- generate.traits(taxa, trait.num, tr, trait.evol.rate)
    
    ### Step 3: Simulate constant rate of preservation
    
    fossils.uni.dupl <- FossilSim::sim.fossils.poisson(rate = fossilisation.rate, tree = tr)
    plot(fossils.uni.dupl, tr, strata = bins, show.strata = TRUE)
    
    ### Step 4: Simulate biogeography on tree
    # assumption: approach assumes migration does not influence tree shape
    # calculate threshold values for number of taxa in each geographic area
    # loop to simulate biogeography as a binary character under the Mk model
    # loop keeps track of number of attempts, exits if iteration.limit is reached
    # inputs tree and migration rate
    # outputs fossil.biogeographic.area, object of type double
    
    #TO DO: integrate resetting iteration.count before 'if'
    
    number_of_tips <- length(traits$sp)
    L <- round(sum(threshold*number_of_tips))
    H <- sum(number_of_tips-L)
    iteration.count <- 0 #always set at 0, resetting it
    #  biogeography.stuck = FALSE
    while(fossils.in.area1 < L || fossils.in.area1 > H) {
      
      ## Running the biogeography simulation
      fossil.biogeographic.area <- FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = migration.rate)
      ## Updating the number of fossils
      fossils.in.area1 <- sum(fossil.biogeographic.area == '1')
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
      if (biogeography.stuck.count >= num.rep*2) {
        stop("Stuck in deep biogeography mud")
      }      
    }
    else
    {
      break
    }
    
  }
  ### Step 5: Simulate biased sampling
  # associate high and low sampling with biogeographical areas in fossil.biogeographic.area [input]
  # simulate biased sampling on tree [input]
  # output is fossils.bio = fossil taxa and respective ages when sampling is biased
  
  rates <- translate.states(fossil.biogeographic.area, low.sampling, high.sampling)
  
  fossils.bio.dupl <- FossilSim::sim.fossils.poisson(rates, tree = tr)
  
  # colours
  fossil.colours = sapply((fossil.biogeographic.area[unlist(sapply(fossils.bio.dupl$sp, function(i) which(taxa$sp == i)))]), function(j) if(j == 1) fossil.colour1 else fossil.colour2)
  plot(fossils.bio.dupl, tr, strata = bins, show.strata = TRUE, fossil.col = fossil.colours)
  
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
  my_geography <- fossil.biogeographic.area
  #my_geography <- data.frame(fossil.biogeographic.area, taxa$sp) ## change to adding bins not taxa numbers
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
