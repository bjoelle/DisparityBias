#plot(simulations[[1]]$matrix[[1]][,1],simulations[[1]]$matrix[[1]][,2])

# x = trait 1, y = trait 2
# color black = all

# color 3 = biased, area 1
# color 4 = biased, area 2

#all = 1
#area 1 = 2
#area 2 = 3


library(ggplot2)
library(plyr)

df = data.frame(trait1 = c(), trait2 = c(), area = c())

# plot all true values
df = rbind(df, data.frame(trait1 = simulations[[1]]$matrix[[1]][,1][simulations[[1]]$subsets$all_species$elements[,1]], 
                          trait2 = simulations[[1]]$matrix[[1]][,2][simulations[[1]]$subsets$all_species$elements[,1]],
                          area = "all species"))
# plot area 1
df = rbind(df, data.frame(trait1 = simulations[[1]]$matrix[[1]][,1][simulations[[1]]$subsets$area_0$elements[,1]], 
                          trait2 = simulations[[1]]$matrix[[1]][,2][simulations[[1]]$subsets$area_0$elements[,1]],
                          area = "all area 0"))

# plot area 2
df = rbind(df, data.frame(trait1 = simulations[[1]]$matrix[[1]][,1][simulations[[1]]$subsets$area_1$elements[,1]], 
                          trait2 = simulations[[1]]$matrix[[1]][,2][simulations[[1]]$subsets$area_1$elements[,1]],
                          area = "all area 1"))

# plot biased towards area 0
df = rbind(df, data.frame(trait1 = simulations[[1]]$matrix[[1]][,1][simulations[[1]]$subsets$bias_0_sample$elements[,1]], 
                          trait2 = simulations[[1]]$matrix[[1]][,2][simulations[[1]]$subsets$bias_0_sample$elements[,1]],
                          area = "biased sampling area 0"))

find_hull <- function(df) df[chull(df$trait1, df$trait2), ]
hulls <- ddply(df, "area", find_hull)

MorphoPlot=ggplot(df, aes(x = trait1, y = trait2, fill=area)) + 
  geom_point(data = subset(df, (area %in% c("all area 1", "all area 0"))), aes(colour = factor(area)), alpha = 0.7, position=position_dodge((width=0.1)))+
  geom_polygon(data = hulls, alpha = 0.2) +
  theme_classic()

pdf(file = paste0(outdir, "MorphoPlot.pdf"), width = wd, height = (ht*3))
par(mfcol=c(1, 3))
print(MorphoPlot)
dev.off()
