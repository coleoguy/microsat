dat <- read.csv("results/busco.results.csv")
hist(dat$score, breaks=25,
     xlab="BUSCO Score",
     main = "Frequency of BUSCO Scores",
     col=rgb(.31,0,0,.5))
