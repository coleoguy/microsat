# read in the results
dat <- read.csv("../results/cent.vs.rate.csv")
dat$diff <- dat$rate.mon-dat$rate.hol
plot(density(dat$diff), xlim=c(-.1,.3),
     cex.lab=.75,main="",
     xlab=expression(paste(sigma^2, " difference (monocentric-holocentric)")))
polygon(density(dat$diff),col=rgb(250,159,181, maxColorValue = 255))
abline(v=0, lty=2,col="gray")
#export pdf at 4.3" x 4.3"


