#load in library for beeswarm
library(beeswarm)

#load in order rate data
order.rates <- read.csv("../results/order.rates.csv")[,-1]

#tidy up the data
colnames(order.rates) <- c("Lepidoptera", "Hymneoptera", "Diptera", 
                           "Coleoptera", "Hemiptera")
results <- data.frame(as.vector(unlist(order.rates)),
                      rep(colnames(order.rates), each = 100))
colnames(results) <- c("rates", "order")

#create the beeswarm plot with the 
beeswarm(log(results$rates) ~ results$order, las = 1,
         xlab = "Order",
         ylab = "Microsatellite Rates (bp change per MY)",
         col = rgb(250, 159, 181, 100,
                   maxColorValue = 255), pch = 16,
         method = "swarm", cex.main = .9,
         cex = .8, spacing = .3, cex.lab = .8, cex.axis = .5
)

#export as pdf 4" x 4"