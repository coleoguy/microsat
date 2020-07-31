#load in data
results <- read.csv("../../results/gsz.content.csv")
str <-read.csv("../../data/str.content.csv")

#put gsz and microsatellite data into different units
str$all/1000000 -> msat

#plot the microsatellite content and genome size
plot(msat~str$gsz,
     ylab = "Microsatellite Content (Mbp)",
     xlab = "Genome Size (Mbp)",
     pch = 16,
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))
lines(x=c(0,2500), lty=2,
      y=c(mean(results$intercept)/1000000,
          (mean(results$beta.content)*2500 + mean(results$intercept))/1000000))
#save as pdf 4.3" x 4.3"
