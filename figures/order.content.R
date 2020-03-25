#set wd as figures
#load in libraries
library(ggplot2)

#read in the microsatellite data
dat.mic <- read.csv("../results/ssr.inference/micRocounter_results_TII.csv", 
                    as.is = T, row.names = 4)

#plot content for orders
ggplot(data = dat.mic, legend = F, aes(order, log(bp.Mbp))) + 
  theme_classic() +
  labs(x = "Order", y = "log(bp.Mbp)") +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width=0.3),
               color = rgb(250, 159, 181, 100, 
                         maxColorValue = 255))+ 
  theme(legend.position = "none") 
  
#export as pdf 8" x 6"