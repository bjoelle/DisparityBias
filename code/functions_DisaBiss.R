
simulation.pipeline <- function(birth, death, tips, trait.num, trait.evol.rate, fossilisation.rate, migration.events, low.sampling, high.sampling, bins, fossil.colour1, fossil.colour2, iteration, variable, variable_i){
  
  pdf(paste0(outdir, "simulated_data_", variable, "_", variable_i, "_", iteration, ".pdf"), height = 11, width = 8.5)
  
  # tree with migration events*
  out = joined_trees(1, tips, migration.events, birth, death)
  tree = out[[1]]
  area = tree$area
  
  # clarify taxonomy
  taxa <- FossilSim::sim.taxonomy(tree = tree, beta = 1)
  taxa$area = 0
  taxa$col = 0
  for(i in c(1:length(taxa$edge))){
    if(taxa$mode[i] == "r" || taxa$mode[i] == "o") { taxa$area[i] = 1; taxa$col[i] = fossil.colour1; next } 
    node = taxa$edge[i]
    taxa$area[i] = area[which(tree$edge[,2] == node)]
    if(taxa$area[i] == 1) 
      taxa$col[i] = fossil.colour1
    else taxa$col[i] = fossil.colour2
  }

  ### Step 2: Simulate "true" disparity
  traits <- generate.traits(taxa, trait.num, tr, trait.evol.rate)
  
  ### Step 3: Simulate constant rate of preservation
  fossils.uni.dupl <- FossilSim::sim.fossils.poisson(rate = fossilisation.rate, taxonomy = taxa)
  plot(fossils.uni.dupl, tree, strata = bins, show.strata = TRUE)
  
  ## Low sampling in area 1, high sampling in area 0 #TODO check this statement is correct
  sampling.rate.0 <- translate.states.0(taxa$area, low.sampling, high.sampling)
  fossils.bias.0.dupl <- FossilSim::sim.fossils.poisson(sampling.rate.0, taxonomy = taxa)
  
  # colourful plots
  fossil.colours.0 <- taxa$col[sapply(fossils.bias.0.dupl$edge, function(i) which(taxa$edge == i))]
  
  plot(fossils.bias.0.dupl, tree, strata = bins, show.strata = TRUE, fossil.col = fossil.colours.0, rho = 0)

  ## Low sampling in area 0, high sampling in area 1
  sampling.rate.1 <- translate.states.1(taxa$area, high.sampling, low.sampling)
  fossils.bias.1.dupl <- FossilSim::sim.fossils.poisson(sampling.rate.1, taxonomy = taxa)
  
  # colourful plots
  fossil.colours.1 <- taxa$col[sapply(fossils.bias.1.dupl$edge, function(i) which(taxa$edge == i))]
  plot(fossils.bias.1.dupl, tree, strata = bins, show.strata = TRUE, fossil.col = fossil.colours.1, rho = 0)
  
  dev.off()
  
  ### Step 6: Bin fossils and match traits with species & bins
  # assumption: no extant samples simulated or sampled, although some fossil species may be extant 
  # calculate bin max/min ages based on tree [input] and number of bins
  
  max.age <- FossilSim::tree.max(tree)
  int.ages <- seq(0, max.age, length = bins + 1)
  
  ###### bin taxa
  # bin all taxa - the following allows us to explore what happens when we have all species trait values for a given interval or area
  bin.all <- bin.taxa(taxa, 3, max.age)
  fossils.all.binned <- FossilSim::sim.interval.ages(bin.all, max.age = max.age, strata = bins, use.species.ages = FALSE)
  
  # bin fossils for unbiased sampling set
  fossils.uni.binned <- FossilSim::sim.interval.ages(fossils.uni.dupl, tree, max.age = max.age, strata = bins, use.species.ages = FALSE)
  # bin fossils for biased sampling set
  fossils.bias.0.binned <- FossilSim::sim.interval.ages(fossils.bias.0.dupl, tree, max.age = max.age, strata = bins, use.species.ages = FALSE)
  fossils.bias.1.binned <- FossilSim::sim.interval.ages(fossils.bias.1.dupl, tree, max.age = max.age, strata = bins, use.species.ages = FALSE)
  
  bias.0 <- int.assign(fossils.bias.0.binned, int.ages)
  bias.1 <- int.assign(fossils.bias.1.binned, int.ages)
  uni <- int.assign(fossils.uni.binned, int.ages)
  # all <- int.assign(fossils.all.binned, int.ages)
  # all$area = sapply(all$sp, function(i) taxa[which(taxa$sp == i),]$area)

  # Grabbing just trait values
  trait.space <- traits[, c("trait1", "trait2")]
  
  uni.sample.int1 <- subset(uni$sp, uni$int == "1")
  uni.sample.int1 <- uni.sample.int1[!duplicated(uni.sample.int1)] #removing duplicates
  uni.sample.int2 <- subset(uni$sp, uni$int == "2")
  uni.sample.int2 <- uni.sample.int2[!duplicated(uni.sample.int2)] #removing duplicates
  # uni.sample.int3 <- subset(uni$sp, uni$int == "3")
  # uni.sample.int3 <- uni.sample.int3[!duplicated(uni.sample.int3)] #removing duplicates
  
  bias.0.sample.int1 <- subset(bias.0$sp, bias.0$int == "1")
  bias.0.sample.int1 <- bias.0.sample.int1[!duplicated(bias.0.sample.int1)]
  bias.0.sample.int2 <- subset(bias.0$sp, bias.0$int == "2")
  bias.0.sample.int2 <- bias.0.sample.int2[!duplicated(bias.0.sample.int2)]
  # bias.0.sample.int3 <- subset(bias.0$sp, bias.0$int == "3")
  # bias.0.sample.int3 <- bias.0.sample.int3[!duplicated(bias.0.sample.int3)]

  bias.1.sample.int1 <- subset(bias.1$sp, bias.1$int == "1")
  bias.1.sample.int1 <- bias.1.sample.int1[!duplicated(bias.1.sample.int1)]
  bias.1.sample.int2 <- subset(bias.1$sp, bias.1$int == "2")
  bias.1.sample.int2 <- bias.1.sample.int2[!duplicated(bias.1.sample.int2)]
  # bias.1.sample.int3 <- subset(bias.1$sp, bias.1$int == "3")
  # bias.1.sample.int3 <- bias.1.sample.int3[!duplicated(bias.1.sample.int3)]
  
  ## Creating the group vector for dispRity
  my.groups <- list(
    # ## All the species
    # "all_species" = subset(all$sp, all$int == "2"),
    # ## All species in location 1
    # #"area_0" = subset(all$sp, all$int == "2")[(subset(all$sp, all$int == "2") %in% which(fossil.biogeographic.area == 0))],
    # "area_0" = subset(all$sp, all$int == "2" & all$area == "1"),
    # ## All species in location 2
    # #"area_1" = subset(all$sp, all$int == "2")[(subset(all$sp, all$int == "2") %in% which(fossil.biogeographic.area == 1))],
    # "area_1" = subset(all$sp, all$int == "2" & all$area == "2"), 
    
    ## The uniform sampled group
    "uni_sample.int1" = uni.sample.int1,
    ## The biased sampled group
    "bias_0_sample.int1" = bias.0.sample.int1,
    ## The biased sampled group
    "bias_1_sample.int1" = bias.1.sample.int1,
    ##int2
    "uni_sample.int2" = uni.sample.int2,
    "bias_0_sample.int2" = bias.0.sample.int2,
    "bias_1_sample.int2" = bias.1.sample.int2
    ##int3
    # "uni_sample.int3" = uni.sample.int3,
    # "bias_0_sample.int3" = bias.0.sample.int3,
    # "bias_1_sample.int3" = bias.1.sample.int3

    # ## unif species sampling
    # "uni_species" = sample(subset(all$sp, all$int == "2"), 20)
  )
  
  ## Creating a dispRity object that contains the trait space and the groups
  disp.groupings <- custom.subsets(data = as.matrix(trait.space),
                                 group = my.groups)
  #TG: ignore the warning (or read it to know what it just did ;) - but nothing bad happening here)
  
  return(disp.groupings)
}

# generate new file for storing traits with taxa in it already [input]
# simulate trait.num number of traits and append to traits file [output]
generate.traits <- function(taxa, trait.num, tr, trait.evol.rate){
  traits <- taxa
  for(i in 1:trait.num){
    tmp <- FossilSim::sim.trait.values(init = 5, taxonomy = taxa, model = "BM", v = trait.evol.rate, min.value = 0)
    traits <- cbind(traits, tmp)
    colnames(traits)[ncol(traits)] <- paste0("trait",i)
  }
  return(traits)
}


# associate high and low sampling with biogeographical areas in fossil.biogeographic.area [input]
translate.states.0 <- function(fossil.biogeographic.area, low.sampling, high.sampling) sapply(fossil.biogeographic.area, function(t) if(t == 1) low.sampling else high.sampling)
translate.states.1 <- function(fossil.biogeographic.area, low.sampling, high.sampling) sapply(fossil.biogeographic.area, function(t) if(t == 0) low.sampling else high.sampling)

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

#### Analysis
disparity.analysis <- function(simulations, analysis = "sum of variances"){
  
  ### Sum of variances
  if(analysis == "sum of variances"){
    # Measure the disparity on the output using lapply (applying a function to a list)
    disparity <- lapply(simulations, dispRity, metric = c(sum, variances))
    title = "Sum of Variances"
  } else if (analysis == "pairwise distance"){
    disparity <- lapply(simulations, dispRity, metric = c(median, pairwise.dist))
    title = "Median pairwise distances"
  } else if (analysis == "centroids") {
    disparity <- lapply(simulations, dispRity, metric = c(median, centroids), centroid = 0)
    title = "Median distance from centroids"
  }else if (analysis == "sum of ranges") {
    disparity <- lapply(simulations, dispRity, metric = c(sum, ranges))
    title = "sum of Ranges"
  }
  
  ####
  # Extract the disparity values (the point estimates explained above)
  point.estimates <- lapply(disparity, get.disparity)
  
  ## Combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
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
  
  # observed_disparity - mean(null_disparity)
  null_disparity <- mean(unlist(results.table[,1]))
  
  # sampling regime
  result$sampling <- as.factor(result$sampling)
  # disparity
  result$values <- as.numeric(result$values) - null_disparity
  
  
  #plots
  p <- ggplot(result, aes(x= factor(sampling, levels= c('bias_0_sample.int1', 'bias_1_sample.int1', 'uni_sample.int1', 
                                                        'bias_0_sample.int2', 'bias_1_sample.int2', 'uni_sample.int2')), y=values)) +
    labs(title = title, x = "sampling") +
    geom_boxplot()

  return(p)
}

perc.intervalle <- function(results.table){
  
  # conversion en matrice numérique
  if(is.list(results.table)) {
    results.table <- matrix(
      unlist(results.table),
      nrow = 10,
      ncol = 6,
      byrow = FALSE
    )
  }
  colnames(results.table) <- c("uni_sample.int1", "bias_0_sample.int1", "bias_1_sample.int1",
                                "uni_sample.int2", "bias_0_sample.int2", "bias_1_sample.int2")
  results.table <- as.data.frame(results.table)
  results.table[results.table == "logical,0"] <- NA
  for(col in colnames(results.table)) {
    results.table[[col]] <- as.numeric(results.table[[col]])
  }
  
  # mean for each bin
  results <- data.frame(
    # int1
    uni_mean_int1 <- mean(results.table$uni_sample.int1, na.rm = TRUE),
    bias0_mean_int1 <- mean(results.table$bias_0_sample.int1, na.rm = TRUE),
    bias1_mean_int1 <- mean(results.table$bias_1_sample.int1, na.rm = TRUE),
    
    # int2
    uni_mean_int2 <- mean(results.table$uni_sample.int2, na.rm = TRUE),
    bias0_mean_int2 <- mean(results.table$bias_0_sample.int2, na.rm = TRUE),
    bias1_mean_int2 <- mean(results.table$bias_1_sample.int2, na.rm = TRUE)
  )
    
    # percentage of disp difference
    results$uni_perc_diff <- ((results$uni_mean_int2 - results$uni_mean_int1) / results$uni_mean_int1) * 100
    results$bias0_perc_diff <- ((results$bias0_mean_int2 - results$bias0_mean_int1) / results$bias0_mean_int1) * 100
    results$bias1_perc_diff <- ((results$bias1_mean_int2 - results$bias1_mean_int1) / results$bias1_mean_int1) * 100
    
    colnames(results) <- c("uni_mean_int1", "bias0_mean_int1", "bias1_mean_int1", "uni_mean_int2", 
                           "bias0_mean_int2", "bias1_mean_int2", "uni_perc_diff", "bias0_perc_diff", "bias1_perc_diff")
  
  return(results)
}



