#Plot for Microsat Content in bp in comparison to BUSCO score

#read in the BUSCO results
dat.busco <- read.csv("results/busco.results.csv", as.is = T)

#read in the micRocounter results
dat.msc.bp <- read.csv("results/micRocounter_results.csv", as.is = T)

#only keep the busco scores that are good quality
good <- dat.busco$good.qual

#plot the total microsatellite content the BUSCO scores for
#only those genomes that are good quality
plot(x = dat.busco$score[good],
     y = jitter(dat.msc.bp$all[good], factor=90),
     pch=16,
     cex=0.4,
     col=rgb(0,.1,1,.4),
     xlab= "BUSCO Score",
     ylab= "Microsat Content (bp)"
     )

#export image as pdf 4" x 6"