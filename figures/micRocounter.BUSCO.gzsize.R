#Plot for Microsat Content in bp in comparison to BUSCO score relative to genome size                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

#read in BUSCO results
dat.busco <- read.csv("../results/busco/busco.results.csv", as.is = T)

#read in the micRocounter results
dat.msc.bp <- read.csv("../results/ssr.inference/micRocounter_results_TII.csv", as.is = T)

#only keep the BUSCO scores that are good quality
good <- dat.busco$good.qual

#plot the total microsatellite content corrected for genome size
#the BUSCO scores for only those genomes that are good quality
plot(x = dat.busco$score[good]*100,
     y = jitter((dat.msc.bp$all[good]/dat.msc.bp$gsz[good])*1000000, 
                factor=90), 
     pch=12,
     cex=0.4,
     col=rgb(250, 159, 181, 100,
             maxColorValue = 255),
     xlab= "BUSCO Score",
     ylab= "Microsat Content Per Mb"
)

#export image as pdf 4.3", 4.3"
