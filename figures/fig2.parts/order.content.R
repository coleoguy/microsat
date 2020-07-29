#load in library for beeswarm
library(beeswarm)

#load in microsatellite data
dat.mic <- read.csv("../../results/micRocounter_results_TII.csv")

orders <- sort(as.character(unique(dat.mic$order)))
nums <- c(3, 12, 78, 1, 15, 68, 22, 2)
labels <- paste(orders, " (", nums, ")", sep="")
#create the beeswarm plot with the
par(mar=c(7,4,4,4))
beeswarm(log(dat.mic$bp.Mbp) ~ dat.mic$order, las = 2,
         xlab = "",
         ylab = "Microsatellite Content (bp/Mbp)",
         col = rgb(250, 159, 181, 100,
                   maxColorValue = 255), pch = 16,
         labels=labels,
         method = "swarm", cex.main = .9,
         cex = .95, spacing = .3, cex.lab = 1, cex.axis = .8
)

#export as pdf 6" x 4"
