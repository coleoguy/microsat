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
par(mar=c(7,4,4,4))
beeswarm(log(results$rates) ~ results$order, las = 2,
         xlab = "",
         ylab = "Microsatellite Rates (bp change per MY)",
         col = rgb(250, 159, 181, 100,
                   maxColorValue = 255), pch = 16,
         method = "swarm", cex.main = .9,
         cex = .95, spacing = .3, cex.lab = 1, cex.axis = 1
)

#export as pdf 6" x 6"
