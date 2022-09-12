#plot(simulations[[1]]$matrix[[1]][,1],simulations[[1]]$matrix[[1]][,2])

# x = trait 1, y = trait 2
# color black = all

# color 3 = biased, area 1
# color 4 = biased, area 2

#all = 1
#area 1 = 2
#area 2 = 3


library(ggplot2)

df = data.frame(trait1 = c(), trait2 = c(), area = c())

# plot all true values
df = rbind(df, data.frame(trait1 = simulations[[1]]$matrix[[1]][,1], 
                          trait2 = simulations[[1]]$matrix[[1]][,2],
                          area = "all fossils"))
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


ggplot(df, aes(x = trait1, y = trait2)) + geom_point(aes(colour = factor(area)))
