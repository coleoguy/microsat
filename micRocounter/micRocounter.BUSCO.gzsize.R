#Plot for Microsat Content in bp in comparison to BUSCO score relative to genome size                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

#read in BUSCO results
dat.busco <- read.csv("results/busco.results.csv", as.is = T)

#read in the micRocounter results
dat.msc.bp <- read.csv("results/micRocounter_results.csv", as.is = T)

#only keep the BUSCO scores that are good quality
good <- dat.busco$good.qual

#plot the total microsatellite content corrected for genome size
#the BUSCO scores for only those genomes that are good quality
plot(x = dat.busco$score[good]*100,
     y = jitter((dat.msc.bp$all[good]/dat.msc.bp$gsz[good])*1000000, 
                factor=90), 
     pch=16,
     cex=0.4,
     col=rgb(0,.1,1,.4),
     xlab= "BUSCO Score",
     ylab= "Microsat Content Per Mb"
)

#export image as pdf 4", 6"
