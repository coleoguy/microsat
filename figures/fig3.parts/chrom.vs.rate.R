#load in data
dat.rates.chrom <- read.csv("../../data/dat.rates.chrom.csv")

#plot the microsatellite content in bp/Mbp and the diploid chromosome number
plot(dat.rates.chrom$rates.evol~dat.rates.chrom$diploid.num,
     xlab = "Diploid Chromosome Number",
     ylab = "Tip Rate (bp change per MY)",
     pch = 16,
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

# export pdf at 4.3" x 4.3"
