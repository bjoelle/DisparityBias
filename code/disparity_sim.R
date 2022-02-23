

set.seed(125)

### Step 1: Simulate tree 

birth <- 0.1
death <- 0.05 
tips <- 100

tr <- TreeSim::sim.bd.taxa(n = tips, 1, birth, death)[[1]]
plot(tr)


### Step 2: Simulate "true" disparity
# current assumption: the trait value for each branch (i.e. each species) is the value at the end of the branch - decision made to maximise diffs between species
# current assumption: bifurcating speciation = each branch is a species

# store trait values with species
taxa <- FossilSim::sim.taxonomy(tr, beta = 1) # how you define morphotaxa with respect to the tree
traits <- taxa

trait_num <- 2 # number of traits we are simulating
v <- 1 # rate of trait evolution

# loop to similate trait_num number of traits
for(i in 1:trait_num){
  tmp <- FossilSim::sim.trait.values(init = 5, tree = tr, model = "BM", v = v, min.value = 0)
  traits <- cbind(traits, tmp)
  colnames(traits)[ncol(traits)] <- paste0("trait",i)
}


## TO DO: add stasis?


### Step 3: Simulate constant rate of preservation

rate <- 0.2 # rate of fossilisation

fossils <- FossilSim::sim.fossils.poisson(rate = rate, tree = tr)

bins <- 3 # number of time bins

plot(fossils, tr, strata = bins, show.strata = TRUE)


### Step 4: Simulate biogeography on tree
# assumption: approach assumes migration does not influence tree shape

rate.bio = 0.04 # migration rate 
traits.bio = FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = rate.bio) # simulating biogeography using Mk model


### Step 5: Simulate biased sampling
low = 0.0015 # low sampling area
high = 0.5 # high sampling area
# associate high and low sampling with biogeographical areas
translate.states = function(traits.bio, low, high) sapply(traits.bio, function(t) if(t == 1) low else high)
rates = translate.states(traits.bio, low, high)
# simulate biased sampling
fossils.bio = FossilSim::sim.fossils.poisson(rates, tree = tr)
plot(fossils.bio, tr, strata = bins, show.strata = TRUE)


### Step 6: Bin fossils and match traits with species & bins
max.age = FossilSim::tree.max(tr)
# define interval ages
int.ages <- seq(0, max.age, length = bins + 1)

# assumption: no extant samples simulated or sampled, although some fossil species may be extant 

# bin fossils for unbiased sampling set
fossils.binned <- FossilSim::sim.interval.ages(fossils, tr, max.age = max.age, strata = bins, use.species.ages = TRUE)
# bin fossils for biased sampling set
fossils.bio.binned <- FossilSim::sim.interval.ages(fossils.bio, tr, max.age = max.age, strata = bins, use.species.ages = TRUE)

# create a new data.frame for disparity analyses based on sampled species in each bin, puts data in format for disparity analysis
disparity.df <- function(traits, fossils, interval.ages){
  if(identical(fossils$hmin, fossils$hmax))
    stop("fossils must be binned!")
  disp <- data.frame()
  for(i in 1:(bins - 1)){
    hmin <- interval.ages[i]
    hmax <- interval.ages[i+1]
    tmp <- subset(fossils, hmin == interval.ages[i])
    tmp <- unique(tmp)
    for(j in tmp$sp){
      tmp2 <- data.frame(bin = i, sp = j, bin.midpoint = mean(c(hmin, hmax)))
      for(k in 1:trait_num){
        t <- traits[which(traits$sp == j),][paste0("trait",k)][[1]]
        tmp2 <- cbind(tmp2, data.frame(tc = t))
        colnames(tmp2)[ncol(tmp2)] <-  paste0("trait",k)
      }
      disp <- rbind(disp, tmp2) 
    }
  }
  disp
}

disp <- disparity.df(traits, fossils.binned, int.ages) ## ? sampled fossils with uniform sampling
disp.bio <- disparity.df(traits, fossils.bio.binned, int.ages) ## ? sampled fossils with biased sampling

#h1 <- hist(traits$trait1)
#h2 <- hist(disp$trait1, breaks = h1$breaks)
#h3 <- hist(disp.bio$trait1, breaks = h1$breaks)

#plot( h1, col=rgb(0,0,1,1/4), xlab = "Trait values", main = NULL)  # first histogram
#plot( h2, col=rgb(1,0,0,1/4), add=T)  # second
#plot( h3, col=rgb(0,0.8,0.6,1/4), add=T)  # third

#legend(x = "topright", legend = c("True disparity", "Uniform sampling","Biased sampling"),
#       fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), rgb(0,.5,.5,1/4)) )

plot( (((traits$start - traits$end) / 2) + traits$start), traits$trait1, col = rgb(0,0,1,1/4), pch = 19)
points(disp$bin.midpoint, disp$trait1, col = rgb(1,0,0,1/4), pch = 19)
points(disp.bio$bin.midpoint, disp.bio$trait1, col = rgb(0,.5,.5,1/4), pch = 19)


### Step 7: Compute disparity metrics
# concatenate matrices
concat_matrix <- matrix(nrow = (nrow(traits)+nrow(disp)+nrow(disp.bio)), ncol = trait_num)
for(i in 1:trait_num){
  t <- paste0("trait",i)
  true <- as.matrix(traits[,t])
  uni <- as.matrix(disp[,t])
  bias <- as.matrix(disp.bio[,t])
  concat_matrix[,i] <- rbind(true, uni, bias)
}

head(concat_matrix)

# keys - simulated vectors may have variable lengths 
true.mx <- length(traits$trait1)
uni.mn <- true.mx + 1
uni.mx <- true.mx + length(disp$trait1)
bias.mn <- uni.mx + 1
bias.mx <- uni.mx + length(disp.bio$trait1)

# ordinating the matrices with traits ### may not need to as traits are independent
ordin.all <- prcomp(concat_matrix)
ordinated_all <- ordin.all$x
rownames(ordinated_all) <- c(1:bias.mx)

# compute metric
library(dispRity)
disparity_data <- dispRity::dispRity.per.group(ordinated_all,
                                     list(trueO = c(1:true.mx), uni = c(uni.mn:uni.mx), bias = c(bias.mn:bias.mx)),
metric = c(median,centroids)) #centroids

disparity_data

plot(disparity_data)




