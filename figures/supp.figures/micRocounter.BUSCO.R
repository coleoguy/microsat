#Plot for Microsat Content in bp in comparison to BUSCO score

#read in the BUSCO results
dat.busco <- read.csv("../../results/busco.results.csv", as.is = T)

#read in the micRocounter results
dat.msc.bp <- read.csv("../../results/micRocounter_results.csv", as.is = T)


#plot the total microsatellite content the BUSCO scores for
#only those genomes that are good quality
plot(x = dat.busco$score,
     y = dat.msc.bp$all,
     pch=16,
     cex=0.4,
     col=rgb(250, 159, 181, 100,
             maxColorValue = 255),
     xlab= "BUSCO Score",
     ylab= "Microsat Content (bp)")
abline(v=.9,col="black",lty=2)

#export image as pdf 6" x 4"


