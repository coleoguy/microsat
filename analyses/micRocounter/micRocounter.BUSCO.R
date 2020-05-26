#Plot for Microsat Content in bp in comparison to BUSCO score

#read in the BUSCO results
dat.busco <- read.csv("../../results/busco/busco.results.csv", as.is = T)

#read in the micRocounter results
dat.msc.bp <- read.csv("../../analyses/micRocounter/micRocounter_results.csv", as.is = T)


#plot the total microsatellite content the BUSCO scores for
#only those genomes that are good quality
plot(x = dat.busco$score,
     y = dat.msc.bp$all,
     pch=16,
     cex=0.4,
     col=rgb(0,.1,1,.4),
     xlab= "BUSCO Score",
     ylab= "Microsat Content (bp)")
abline(v=.9,col="red",lty=2)

#export image as pdf 4" x 6"



hist(dat.busco$score[!dat.msc.bp$good.qual], breaks=100)
hist(dat.busco$score[dat.msc.bp$good.qual], breaks=100)

