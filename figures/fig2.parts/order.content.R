#load in library for beeswarm
library(beeswarm)

#load in microsatellite data
dat.mic <- read.csv("../results/micRocounter_results_TII.csv")


#create the beeswarm plot with the 
beeswarm(log(dat.mic$bp.Mbp) ~ dat.mic$order, las = 1,
         xlab = "Order",
         ylab = "Microsatellite Content (bp/Mbp)",
         col = rgb(250, 159, 181, 100,
                   maxColorValue = 255), pch = 16,
         method = "swarm", cex.main = .9,
         cex = .8, spacing = .3, cex.lab = .8, cex.axis = .5
)

#export as pdf 6" x 4"
