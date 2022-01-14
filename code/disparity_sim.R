

set.seed(125)

# simulate tree

birth <- 0.1
death <- 0.05 
tips <- 100

tr <- TreeSim::sim.bd.taxa(n = tips, 1, birth, death)[[1]]
plot(tr)

# simulate "true" disparity
# current assumption: the trait value for each branch (i.e. each species) is the value at the end of the branch - decision made to maximise diffs between species
# current assumption: bifurcating speciation = each branch is a species

# store trait values with species
taxa <- FossilSim::sim.taxonomy(tr, beta = 1)
traits <- taxa

trait_num <- 2
v <- 1 # rate of trait evolution

for(i in 1:trait_num){
  tmp <- FossilSim::sim.trait.values(init = 5, tree = tr, model = "BM", v = v, min.value = 0) #removed min.value = 0
  traits <- cbind(traits, tmp)
  colnames(traits)[ncol(traits)] <- paste0("trait",i)
}



## TO DO: add stasis?

# simulate constant rate of preservation

bins <- 3

rate <- 0.2 ### decreased slightly from 0.3 to reduce number of fossils overall
fossils <- FossilSim::sim.fossils.poisson(rate = rate, tree = tr)

plot(fossils, tr, strata = bins, show.strata = TRUE)

# simulate biogeography 
# note approach assumes migration does not influence tree shape

rate.bio = 0.0075 # migration rate ### increased slightly from 0.005, few taxa from area 2 otherwise
traits.bio = FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = rate.bio)

low = 0.015
high = 0.3 #### decreased from 0.5, rate of preservation was super high
translate.states = function(traits.bio, low, high) sapply(traits.bio, function(t) if(t == 0) low else high)
############# are we both making migration into area 2 unlikely, and under-sampling area 2?
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
  #disp <- data.frame(bin = c(), sp = c(), traits = c(), bin.midpoint = c())
  disp <- data.frame(bin = c(), sp = c(), bin.midpoint = c())  
  for(i in 1:(bins - 1)){
    hmin <- interval.ages[i]
    hmax <- interval.ages[i+1]
    tmp <- subset(fossils, hmin == interval.ages[i])
    tmp <- unique(tmp)
    disp <- data.frame()
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

h1 <- hist(traits$trait1)
h2 <- hist(disp$trait1, breaks = h1$breaks)
h3 <- hist(disp.bio$trait1, breaks = h1$breaks)

plot( h1, col=rgb(0,0,1,1/4), xlab = "Trait values", main = NULL)  # first histogram
plot( h2, col=rgb(1,0,0,1/4), add=T)  # second
plot( h3, col=rgb(0,0.8,0.6,1/4), add=T)  # third

legend(x = "topright", legend = c("True disparity", "Uniform sampling","Biased sampling"),
       fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), rgb(0,.5,.5,1/4)) )

plot( (((traits$start - traits$end) / 2) + traits$start), traits$trait1, col = rgb(0,0,1,1/4), pch = 19)
points(disp$bin.midpoint, disp$trait1, col = rgb(1,0,0,1/4), pch = 19)
points(disp.bio$bin.midpoint, disp.bio$trait1, col = rgb(0,.5,.5,1/4), pch = 19)

#Catherine does her whole disparity thing :)
# Catherine *tries* to do the disparity thing

## If I understand correctly, traits$traits is values to calculate true disparity, 
## disp$traits for uniform sampling, disp.bio$traits for biased sampling

# concatenate matrices
concat_matrix <- matrix(nrow = (nrow(traits)+nrow(disp)+nrow(disp.bio)), ncol = trait_num)
for(i in 1:trait_num){
  t <- paste0("trait",i)
  true <- as.matrix(traits[,t])
  uni <- as.matrix(disp[,t])
  bias <- as.matrix(disp.bio[,t])
  concat_matrix[,i] <- rbind(true, uni, bias)#cbind(concat_matrix, tmp)
}

concat_matrix

# keys - simulated vectors may have variable lengths
true.mx <- length(traits$trait1)
uni.mn <- true.mx + 1
uni.mx <- true.mx + length(disp$trait1)
bias.mn <- uni.mx + 1
bias.mx <- uni.mx + length(disp.bio$trait1)

#ordinating the matrices with traits
ordin.all <- prcomp(concat_matrix)
ordinated_all <- ordin.all$x

# separating ordinated matrix ###not necessary
#trueO <- ordinated_all[1:369]
#uniO <- ordinated_all[370:386]
#biasO <- ordinated_all[387:409]

#groups <- list("trueO" = c(1:369), "uni" = c(370:386), "bias" = c(387:409))
#rownames(ordinated_all) <- c(1:409)
rownames(ordinated_all) <- c(1:bias.mx)

#ordin2 <- cbind(ordinated_all, replicate(2,ordinated_all[,1]))

disparity_data <- dispRity::dispRity.per.group(ordinated_all,
                                     list(trueO = c(1:true.mx), uni = c(uni.mn:uni.mx), bias = c(bias.mn:bias.mx)),
metric = c(median,centroids)) #centroids

disparity_data




############# stuff below is from when I ran dispRity previously with my extant datasets
timebin <- list("Ypresian" = c(1:369), "Lutetian" = c(370:386), "Bartonian" = c(387:408))

PCdataCust <- custom.subsets(as.matrix(ordinated_all), group = timebin)
###### Warning: In custom.subsets(PCdata1, group = timebin) :Rownames generated for PCdata1 as seq(1:148)


##### Bootstrap
PCdataBoots <- boot.matrix(PCdataCust, bootstraps = 100, boot.type = "full", rarefaction = TRUE)
# which level do I want to rarefy to?
# how many times to bootstrap?


#### Sum of variances
SumV <- dispRity(PCdataBoots, metric = c(sum, variances))
summary(SumV)

plot(SumV,type = "box")

