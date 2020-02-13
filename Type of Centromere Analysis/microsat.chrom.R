#read in microsatellite data
dat.mic <- read.csv("results/micro.vs.chrom.csv", as.is=T)

#read in chromosome data
dat.chroms <- read.csv("data/chroms/data.invert.csv", as.is=T)

#vector of names from chromosome data
chrom.names <- paste(dat.chroms$Genus, dat.chroms$species)

#loop that finds species diploid chromosome number from our data frame
for(i in 1:nrow(dat.mic)){
  # if species in microcounter matches one in chromosome data 
  if(dat.mic$species[i] %in% chrom.names){
    #store the name in vector hit
    hit <- which(chrom.names == dat.mic$species[i])[1]
    #fill in diploid number for those species that have a match in the 
    #chromosome data
    dat.mic$diploid.num[i] <- dat.chroms$Chromosome.number..female..2N[hit]
  }
}

#make a new data frame that will contain only those species with both 
#microsatellite and chromosome number data
dat.intersect <- dat.mic[complete.cases(dat.mic),]

#store the bp/Mbp microsatellite content
bpMbp <- dat.intersect$all/(dat.intersect$gsz/1000000)

#plot the microsatellite content in bp/Mbp and the diploid chromosome number
plot(bpMbp~dat.intersect$diploid.num,
     xlab = "Diploid Number",
     ylab = "Microsatellite Content (bp/Mbp)")

#save plot as pdf 4" x 6" 

#write the CSV with microsatellite, genome, and diploid chromosome number data
write.csv(dat.intersect, "dat.intersect.csv")
