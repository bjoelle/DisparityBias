

set.seed(125)


### Setting up variables
birth <- 0.1 #birth rate
death <- 0.05 # death rate
tips <- 100 # number of tips in tree
trait_num <- 2 # number of traits we are simulating
v <- 0.5 # rate of trait evolution
rate <- 0.2 # rate of fossilisation
bins <- 3 # number of time bins
rate.bio = 0.02 # migration rate 
fossils_in_area1 <- 0 # setting up parameter for checking spatial split
threshold <- 0.45 # threshold for spatial split between areas 0 and 1
low = 0.0015 # sampling rate for fossils in low sampling area
high = 0.5 # sampling rate for fossils in high sampling area


### Step 1: Simulate tree 
tr <- TreeSim::sim.bd.taxa(n = tips, 1, birth, death)[[1]]
plot(tr)


### Step 2: Simulate "true" disparity
# current assumption: the trait value for each branch (i.e. each species) is the value at the end of the branch - decision made to maximise diffs between species
# current assumption: bifurcating speciation = each branch is a species

# store trait values with species
taxa <- FossilSim::sim.taxonomy(tr, beta = 1) # how you define morphotaxa with respect to the tree
traits <- taxa

# loop to simulate trait_num number of traits
for(i in 1:trait_num){
  tmp <- FossilSim::sim.trait.values(init = 5, tree = tr, model = "BM", v = v, min.value = 0)
  traits <- cbind(traits, tmp)
  colnames(traits)[ncol(traits)] <- paste0("trait",i)
}
check <- FossilSim::as.fossils(taxa, from.taxonomy = FALSE)

## TO DO: add stasis?


### Step 3: Simulate constant rate of preservation

fossils.uni.dupl <- FossilSim::sim.fossils.poisson(rate = rate, tree = tr)
plot(fossils.uni.dupl, tr, strata = bins, show.strata = TRUE)
fossils.uni <- dplyr::distinct(fossils.uni.dupl, sp, .keep_all = TRUE)
plot(fossils.uni, tr, strata = bins, show.strata = TRUE)

### Step 4: Simulate biogeography on tree
# assumption: approach assumes migration does not influence tree shape

#traits.bio = FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = rate.bio) # simulating biogeography using Mk model

# Making sure the number of taxa in each area is nearly equal

## Setting up the parameter of choice

## Lowest fraction of taxa in area 1 permitted
number_of_tips <- length(tmp)
L <- round(sum(threshold*number_of_tips))
H <- sum(number_of_tips-L)

## Run the while loop to get a set of around 100 fossils (+/1 20)
while(fossils_in_area1 < L || fossils_in_area1 > H) {
  ## Running the biogeography simulation
  traits.bio = FossilSim::sim.trait.values(1, tree = tr, model = "Mk", v = rate.bio) # simulating biogeography using Mk model
  ## Updating the number of fossils
  fossils_in_area1 <- sum(traits.bio == '1')
}

### Step 5: Simulate biased sampling
# associate high and low sampling with biogeographical areas
translate.states = function(traits.bio, low, high) sapply(traits.bio, function(t) if(t == 1) low else high)
rates = translate.states(traits.bio, low, high)
# simulate biased sampling
fossils.bio.dupl = FossilSim::sim.fossils.poisson(rates, tree = tr)
plot(fossils.bio.dupl, tr, strata = bins, show.strata = TRUE)
#fossils.bio <- dplyr::distinct(fossils.bio.dupl, sp, .keep_all = TRUE) #removing duplicates
#plot(fossils.bio, tr, strata = bins, show.strata = TRUE)

### Step 6: Bin fossils and match traits with species & bins
max.age = FossilSim::tree.max(tr)
# define interval ages
int.ages <- seq(0, max.age, length = bins + 1)

# assumption: no extant samples simulated or sampled, although some fossil species may be extant 

# bin fossils for unbiased sampling set
fossils.binned <- FossilSim::sim.interval.ages(fossils.uni.dupl, tr, max.age = max.age, strata = bins, use.species.ages = FALSE)
# bin fossils for biased sampling set
fossils.bio.binned <- FossilSim::sim.interval.ages(fossils.bio.dupl, tr, max.age = max.age, strata = bins, use.species.ages = FALSE)

int.assign <- function(fossils, ints){
  if(identical(fossils$hmin, fossils$hmax))
    stop("fossils must be binned!")
  fossils$int <- NA
  for(i in 1:(length(ints) - 1)){
    fossils[which(fossils$hmin == ints[i]),]$int = i
  }
  fossils
}

test = int.assign(fossils.bio.binned, int.ages)

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


### Step 7: Compute disparity metrics
# concatenate matrices
#concat_matrix <- matrix(nrow = (nrow(traits)+nrow(disp)+nrow(disp.bio)), ncol = trait_num)
#for(i in 1:trait_num){
#  t <- paste0("trait",i)
#  true <- as.matrix(traits[,t])
#  uni <- as.matrix(disp[,t])
#  bias <- as.matrix(disp.bio[,t])
#  concat_matrix[,i] <- rbind(true, uni, bias)
#}

#head(concat_matrix)

# keys - simulated vectors may have variable lengths 
#true.mx <- length(traits$trait1)
#uni.mn <- true.mx + 1
#uni.mx <- true.mx + length(disp$trait1)
#bias.mn <- uni.mx + 1
#bias.mx <- uni.mx + length(disp.bio$trait1)

# ordinating the matrices with traits ### may not need to as traits are independent #TG: probably not needed indeed
#ordin.all <- prcomp(concat_matrix) #TG: Note also that here some species are duplicated (the traitspace contains 413 species/samples but this trait space has 443 - the bias/uni sampled ones). This creates a bias in the PCA.
#ordinated_all <- ordin.all$x
#rownames(ordinated_all) <- c(1:bias.mx)

# compute metric - centroids
library(dispRity)
disparity_centr <- dispRity::dispRity.per.group(ordinated_all,
                                     list(trueO = c(1:true.mx), uni = c(uni.mn:uni.mx), bias = c(bias.mn:bias.mx)),
  metric = centroids) #variance
# metric = variances)
disparity_centr

plot(disparity_centr)

#compute metric - sum of variances
####c(mean,variance))
disparity_var <- dispRity::dispRity.per.group(ordinated_all,
                                               list(trueO = c(1:true.mx), uni = c(uni.mn:uni.mx), bias = c(bias.mn:bias.mx)),
                                               metric = variances) #variance
# metric = variances)
disparity_var

plot(disparity_var)

plot(disparity_centr, type = "preview")









## Example script using a different dispRity pipeline

## Rename some variables
my_geography <- traits.bio
#my_geography <- data.frame(traits.bio, taxa$sp) ## change to adding bins not taxa numbers
my_trait_space <- traits[, c("trait1", "trait2")]

## Creating the group vector for dispRity
my_groups <- list(## All the species
                  "all_species" = 1:nrow(my_trait_space),
                  ## All species in location 1
                  "area_0" = which(my_geography == 0),
                  ## All species in location 2
                  #"area_1" = which((my_geography$traits.bio == 1) & (my_geography$taxa.sp == 2)), #change column name
                  ## The uniform sampled group
                  "uni_sample" = subset(disp$sp, disp$bin == "2"), #TG: assuming that that $sp column contains the species ID/row number in traits[, c("trait1", "trait2")]
                  ## The biased sampled group
                  "bias_sample" = subset(disp.bio$sp, disp.bio$bin == "2"))

## Creating a dispRity object that contains the trait space and the groups
my_groupings <- custom.subsets(data = my_trait_space,
                               group = my_groups)
#TG: ignore the warning (or read it to know what it just did ;) - but nothing bad happening here)

## Calculating disparity on these groups
disparity_sum_var <- dispRity(my_groupings, metric = c(sum, variances))

## Hop
plot(disparity_sum_var)
#TG: note that these are now point estimates (one disparity value per group) hence the absence of variance. You'll get the variance from replicating the simulations (see pseudo code below).
#TG: to just get the values you can use the function get.disparity (see example in the pseudo code below)
get.disparity(disparity_sum_var)



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
my_5_simulations <- repliate(5, simulation.pipeline(birth = X, death = Y, tips = Z, <other_magic_numbers>))

#TG: You can then measure the disparity on the output using lapply (applying a function to a list)
my_sum_variances <- lapply(my_5_simulations, dispRity, metric = c(sum, variances))

#TG: You can then extract the disparity values (the point estimates explained above)
my_point_estimates <- lapply(my_sum_variances, get.disparity)

##TG: And you can combine that into a more reader friendly format (a table!) using the rbind function (bind in rows) applied to this list of lists using do.call
my_results_table <- do.call(rbind, my_point_estimates)
boxplot(my_results_table)
