#set working directory as analyses/order.rates
#read in order rates data
order.rates <- read.csv("../../results/order.rates/order.rates.csv")[,-1]
#adjust column names
colnames(order.rates) <- c("Lepidoptera", "Hymneoptera", "Diptera", "Coleoptera", "Hemiptera")

#plot the order data (boxplot)
boxplot(order.rates, 
        outpch = NA,
        ylim = c(0,0.42),
        xlab = "Order",
        ylab = "Rate Estimates")
stripchart(order.rates, 
           data = order.rates,
           vertical = TRUE, 
           method = "jitter", 
           add = TRUE, 
           pch = 20, 
           col = rgb(250, 159, 181, 100,
                     maxColorValue = 255))

#export as pdf 8" x 6"
