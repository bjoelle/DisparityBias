
####Function 1: Makes x amount of trees with simulated traits####
Tree.Taxa = function(birth, death, tips, trait.num, trait.evol.rate, fossilisation.rate, migration.events, low.sampling, high.sampling, 
                     bins, iteration, variable, variable_i){
  
  # generate tree with migration events
  out = joined_trees(1, tips, migration.events, birth, death)
  tree = out[[1]]
  
  # clarify taxonomy (I don't know what they meant with 'Clarify' -J)
  taxa <- FossilSim::sim.taxonomy(tree = tree, beta = 1)
  taxa$area = 0
  taxa$col = 0
  for(i in c(1:length(taxa$edge))){
    if(taxa$mode[i] == "r" || taxa$mode[i] == "o") { taxa$area[i] = 1; taxa$col[i] = fossil.colour1; next } 
    node = taxa$edge[i]
    taxa$area[i] = tree$area[which(tree$edge[,2] == node)]
    if(taxa$area[i] == 1) 
      taxa$col[i] = fossil.colour1
    else taxa$col[i] = fossil.colour2
  }

  ### Simulate "true" disparity
  traits <- generate.traits(taxa, trait.num, trait.evol.rate)
  
  ### Simulate constant rate of preservation
  fossils.uni.dupl <- FossilSim::sim.fossils.poisson(rate = fossilisation.rate, taxonomy = taxa)

  ## Low sampling in area 1, high sampling in area 2
  sampling.rate.2 <- translate.states(taxa$area, low.sampling, high.sampling, 1)
  fossils.bias.2.dupl <- FossilSim::sim.fossils.poisson(sampling.rate.2, taxonomy = taxa)

  ## Low sampling in area 2, high sampling in area 1
  sampling.rate.1 <- translate.states(taxa$area, low.sampling, high.sampling, 2)
  fossils.bias.1.dupl <- FossilSim::sim.fossils.poisson(sampling.rate.1, taxonomy = taxa)
  ## Save output as RData file.
  Tree.Taxa.Output=list(tree, taxa, traits, fossils.uni.dupl, fossils.bias.2.dupl, fossils.bias.1.dupl)
  save(Tree.Taxa.Output,file=paste0(outdir, "TreeTaxa_", variable, "_",variable_i, "_", iteration, ".RData"))
}


####Function 2: Makes plots of the trees made with function 1####  
Phylogeny.Plots=function(iteration, variable, variable_i, fossil.colour1, fossil.colour2){

  #Load in required RData file
  load(file = paste0(outdir, "TreeTaxa_", variable, "_",variable_i, "_", iteration, ".RData"))
  
  #grab required lists from RData file
  tree=Tree.Taxa.Output[[1]]
  taxa=Tree.Taxa.Output[[2]]
  fossils.uni.dupl=Tree.Taxa.Output[[4]]
  fossils.bias.2.dupl=Tree.Taxa.Output[[5]]
  fossils.bias.1.dupl=Tree.Taxa.Output[[6]]
  
  #Begin pdf process:
  pdf(paste0(outdir, "simulated_data_", variable, "_", variable_i, "_", iteration, ".pdf"), height = 11, width = 8.5)
 
  # allow dots to be coloured on plots:
  fossil.colours.1 <- taxa$col[sapply(fossils.bias.1.dupl$edge, function(i) which(taxa$edge == i))]
  fossil.colours.2 <- taxa$col[sapply(fossils.bias.2.dupl$edge, function(i) which(taxa$edge == i))] 
  #The plots:
  plot(fossils.uni.dupl, tree, strata = bins, show.strata = TRUE)
  plot(fossils.bias.2.dupl, tree, strata = bins, show.strata = TRUE, fossil.col = fossil.colours.2, rho = 0)
  plot(fossils.bias.1.dupl, tree, strata = bins, show.strata = TRUE, fossil.col = fossil.colours.1, rho = 0)

  
  dev.off()
}

  
####Function 3: takes the output from function 1 and make into input fro DispRity ####
DispRity_Input=function(iteration, variable, variable_i){  
  ### Bin fossils and match traits with species & bins
  # assumption: no extant samples simulated or sampled, although some fossil species may be extant 
  # calculate bin max/min ages based on tree [input] and number of bins
  
  #Load in required Rdata file and grab needed lists
  load(file = paste0(outdir, "TreeTaxa_", variable, "_",variable_i, "_", iteration, ".RData"))
  tree=Tree.Taxa.Output[[1]]
  taxa=Tree.Taxa.Output[[2]]
  traits=Tree.Taxa.Output[[3]]
  fossils.uni.dupl=Tree.Taxa.Output[[4]]
  fossils.bias.2.dupl=Tree.Taxa.Output[[5]]
  fossils.bias.1.dupl=Tree.Taxa.Output[[6]]
  
  
  
  max.age <- FossilSim::tree.max(tree)
  int.ages <- seq(0, max.age, length = bins + 1)
  
  ###### bin taxa
  # bin all taxa - the following allows us to explore what happens when we have all species trait values for a given interval or area
  bin.all <- bin.taxa(taxa, 3, max.age) # 3=magic number for nbins, should be defined somewhere
  fossils.all.binned <- FossilSim::sim.interval.ages(bin.all, max.age = max.age, strata = bins, use.species.ages = FALSE)
  
  # bin fossils for unbiased sampling set
  fossils.uni.binned <- FossilSim::sim.interval.ages(fossils.uni.dupl, tree, max.age = max.age, strata = bins, use.species.ages = FALSE)
  # bin fossils for biased sampling set
  fossils.bias.2.binned <- FossilSim::sim.interval.ages(fossils.bias.2.dupl, tree, max.age = max.age, strata = bins, use.species.ages = FALSE)
  fossils.bias.1.binned <- FossilSim::sim.interval.ages(fossils.bias.1.dupl, tree, max.age = max.age, strata = bins, use.species.ages = FALSE)
  
  #NOTE: Bias.0= bias.2. Haven't changed it in this part yet.
  bias.0 <- int.assign(fossils.bias.2.binned, int.ages)
  bias.1 <- int.assign(fossils.bias.1.binned, int.ages)
  uni <- int.assign(fossils.uni.binned, int.ages)
   all <- int.assign(fossils.all.binned, int.ages)
   all$area = sapply(all$sp, function(i) taxa[which(taxa$sp == i),]$area)

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
     "all_species" = subset(all$sp, all$int == "1"),
    # ## All species in location 1
    # #"area_0" = subset(all$sp, all$int == "2")[(subset(all$sp, all$int == "2") %in% which(fossil.biogeographic.area == 0))],
     "area_0" = subset(all$sp, all$int == "1" & all$area == "1"),
    # ## All species in location 2
    # #"area_1" = subset(all$sp, all$int == "2")[(subset(all$sp, all$int == "2") %in% which(fossil.biogeographic.area == 1))],
     "area_1" = subset(all$sp, all$int == "1" & all$area == "2"), 
    
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

####Function 4: Generate trait values. Used in Function 1####
#Has to be adapted to make sense for trilobites

generate.traits <- function(taxa, trait.num, trait.evol.rate){
  traits <- taxa
  for(i in 1:trait.num){
    tmp <- FossilSim::sim.trait.values(init = 5, taxonomy = taxa, model = "BM", v = trait.evol.rate, min.value = 0)
    traits <- cbind(traits, tmp)
    colnames(traits)[ncol(traits)] <- paste0("trait",i)
  }
  return(traits)
}
#### Function 5: Ads low sampling to one subset and high to other? Used in Function 1#### 
# associate high and low sampling with biogeographical areas in fossil.biogeographic.area [input]
translate.states <- function(Area.num, Samp.num, Samp.notNum, Number) sapply(Area.num, function(t) if(t == Number) Samp.num else Samp.notNum)


#### Function 6: turns all taxa into a format useable by FossilSim so that time binning can occur. Used in Function 1####
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

#### Function 7: function to turn sim.interval.ages into defined/numbered time bins. Used in Function 1####
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

#### Function 8: Analysis ####
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
  
  #
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


####Function 9: French?####
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


####Function 10: Analysis of output####
Analysis_Output=function(variable, variable_i){  

  load(file = paste0(outdir, "data_", var, "_", i, "_", ".RData"))
  if(analysis){
    ### Disparity Analysis - these functions return plots
    sumv <- disparity.analysis(simulations, analysis = "sum of variances")
    mpd <- disparity.analysis(simulations, analysis = "pairwise distance")
    mcd <- disparity.analysis(simulations, analysis = "centroids")
    sumr <- disparity.analysis(simulations, analysis = "sum of ranges")
    
    assign(paste0("sumv_", var, "_", i), sumv)
    assign(paste0("mpd_", var, "_", i), mpd)
    assign(paste0("mcd_", var, "_", i), mcd)
    assign(paste0("sumr_", var, "_", i), sumr)
    
    perc_sumv <- perc.intervalle(sumv$plot_env$results.table)
    perc_mpd <- perc.intervalle(mpd$plot_env$results.table)
    perc_mcd <- perc.intervalle(mcd$plot_env$results.table)
    perc_sumr <- perc.intervalle(sumr$plot_env$results.table)
    
    assign(paste0("perc_sumv_", var, "_", i), perc_sumv)
    assign(paste0("perc_mpd_", var, "_", i), perc_mpd)
    assign(paste0("perc_mcd_", var, "_", i), perc_mcd)
    assign(paste0("perc_sumr_", var, "_", i), perc_sumr)
  }
}
