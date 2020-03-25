#set wd as figures
#load in libraries
library(ggplot2)

#read in the microsatellite data
dat.mic <- read.csv("../results/ssr.inference/micRocounter_results_TII.csv",
                    as.is = T, row.names = 4)
true.ticks <- 4:9
print.ticks <- round(exp(true.ticks))
#plot content for orders
ggplot(data = dat.mic, legend = F, aes(order, log(bp.Mbp))) +
  theme_bw() +
  labs(x = "", y = "bp/Mbp") +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width=0.3),
               color = rgb(250, 159, 181, 100,
                         maxColorValue = 255),
             size=.7)+
  scale_y_continuous(breaks=true.ticks, labels=print.ticks) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# export 3.5x3.5
