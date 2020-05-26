#set working directory as analyses/order.rates
#read in order rates data
order.rates <- read.csv("../results/order.rates/order.rates.csv")[,-1]
#adjust column names
colnames(order.rates) <- c("Lepidoptera", "Hymneoptera", "Diptera", "Coleoptera", "Hemiptera")


results <- data.frame(as.vector(unlist(order.rates)),
                     rep(colnames(order.rates), each=100))
colnames(results) <- c("rates", "order")

true.ticks <- c(-7.5,-5,-2.5,0)
print.ticks <- round(exp(true.ticks),digits=3)
#plot content for orders
ggplot(data = results, legend = F, aes(order, log(rates))) +
  theme_bw() +
  #ylim(0,.6) +
  labs(x = "", y = "Rate estimate") +
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
