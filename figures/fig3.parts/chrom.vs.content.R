#read in data
dat.intersect <- read.csv("../../data/dat.intersect.csv")

#store the bp/Mbp microsatellite content
bpMbp <- dat.intersect$all/(dat.intersect$gsz/1000000)

#plot the microsatellite content in bp/Mbp and the diploid chromosome number
plot(bpMbp~dat.intersect$diploid.num,
     xlab = "Diploid Chromosome Number",
     ylab = "Microsatellite Content (bp/Mbp)",
     pch = 16, 
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

#export pdf at 4.3" x 4.3"