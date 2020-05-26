#load in the data
results <- read.csv("../results/gsz.rates.csv")
str <-read.csv("../data/str")

#plot the microsatellite evolution rates and genome size
plot(str$rates ~ str$gsz,
     xlab = "Genome Size (Mbp)",
     ylab = "Tip Rate (bp change per MY)",
     pch = 16,
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

lines(x=c(0,2500), lty=2,
      y=c(mean(results$intercept),
          (mean(results$beta.rates)*2500 +
             mean(results$intercept))))

#export as pdf 4.3" x 4.3"
