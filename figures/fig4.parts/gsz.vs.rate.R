#load in the data
results <- read.csv("../results/gsz.rates.csv")
str <-read.csv("../data/str.rates")

#plot the microsatellite evolution rates and genome size

gsz <- str$gsz/1000000
plot(str$rates ~ gsz,
     xlab = "Genome Size (Mbp)",
     ylab = "Tip Rate (bp change per MY)",
     pch = 16,
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

lines(x=c(0,2500), lty=2,
      y=c(mean(results$intercept)/1000000,
          (mean(results$beta.rates)*2500 + mean(results$intercept))/1000000))

#export as pdf 4.3" x 4.3"
