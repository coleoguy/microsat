#read in the parsed BUSCO results
dat <- read.csv("results/busco.results.csv")

#make a histogram of the BUSCO scores
hist(dat$score, breaks=25,
     xlab="BUSCO Score",
     col=rgb(.31,0,0,.5),
     main = "Disribution of BUSCO Scores")
