#set figures as working directory
# read in the results
dat <- read.csv("../../results/cent.vs.rate.csv")

#calculates and adds column for the difference in monocentric and holocentric rates
dat$diff <- dat$rate.mon-dat$rate.hol

#density plot of the difference in rates
plot(density(dat$diff), xlim = c(-.1,.3),
     cex.lab = .75, main = "",
     xlab = expression(paste(sigma^2, " difference (monocentric-holocentric)")))

#fills in under the curve
polygon(density(dat$diff), col = rgb(250,159,181, maxColorValue = 255))

#places a line at 0 
abline(v = 0, lty = 2,col = "gray")

#export pdf at 4.3" x 4.3"


