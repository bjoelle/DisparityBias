

set.seed(123)

# simulate tree

birth <- 0.1
death <- 0.05 
tips <- 50

tr <- TreeSim::sim.bd.taxa(n = tips, 1, birth, death)[[1]]
plot(tr)

# simulate "true" disparity
# current assumption: the trait value for each branch (i.e. each species) is the value at the end of the branch - decision made to maximise diffs between species
# current assumption: bifurcating speciation = each branch is a species

v = 0.5 # rate of trait evolution
traits <- FossilSim::sim.trait.values(init = 5, tree = tr, model = "BM", v = v, min.value = 0)

# store trait values with species
taxa <- FossilSim::sim.taxonomy(tr, beta = 1)
traits <- cbind(taxa, traits)

## TO DO: add stasis?

# simulate constant rate preservation

bins <- 10

rate <- 0.4
fossils <- FossilSim::sim.fossils.poisson(rate = rate, tree = tr)

plot(fossils, tr, strata = bins, show.strata = TRUE)

# simulate biogeography 
# note approach assumes migration does not influence tree shape

rate.bio = 0.02 # migration rate
traits.bio = FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = rate.bio)

low = 0.05
high = 0.95
translate.states = function(traits.bio, low, high) sapply(traits.bio, function(t) if(t == 0) low else high)

rates = translate.states(traits.bio, low, high)
fossils.bio = FossilSim::sim.fossils.poisson(rates, tree = tr)
plot(fossils.bio, tr, strata = bins, show.strata = TRUE)

### bin fossils and match traits with species & bins
max.age = FossilSim::tree.max(tr)

## no extant samples simulated or sampled, although some fossil species may be extant 
fossils.binned <- FossilSim::sim.interval.ages(fossils, tr, max.age = max.age, strata = bins, use.species.ages = TRUE)

fossils.bio.binned <- FossilSim::sim.interval.ages(fossils.bio, tr, max.age = max.age, strata = bins, use.species.ages = TRUE)

# define interval ages
int.ages <- seq(0, max.age, length = bins + 1)

# create a new data.frame for disparity analyses based on sampled species in each bin

disparity.df <- function(traits, fossils, interval.ages){
  if(identical(fossils$hmin, fossils$hmax))
    stop("fossils must be binned!")
  disp <- data.frame(bin = c(), sp = c(), traits = c())  
  for(i in 1:(bins - 1)){
    hmin <- interval.ages[i]
    tmp <- subset(fossils, hmin == interval.ages[i])
    tmp <- unique(tmp)
    for(j in tmp$sp){
      t <- traits[which(traits$sp == j),]$traits[1]
      disp <- rbind(disp, data.frame(bin = i, sp = j, traits = t)) 
    }
  }
  disp
}

disp <- disparity.df(traits, fossils.binned, int.ages) 
disp.bio <- disparity.df(traits, fossils.bio.binned, int.ages) 

h1 <- hist(traits$traits)
h2 <- hist(disp$traits, breaks = h1$breaks)
h3 <- hist(disp.bio$traits, breaks = h1$breaks)

plot( h1, col=rgb(0,0,1,1/4), xlab = "Trait values", main = NULL)  # first histogram
plot( h2, col=rgb(1,0,0,1/4), add=T)  # second
plot( h3, col=rgb(0,0.8,0.6,1/4), add=T)  # third

legend(x = "topright", legend = c("True disparity", "Uniform sampling","Biased sampling"),
       fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4),rgb(0,.5,.5,1/4)) )


#Catherine does her whole disparity thing :)
