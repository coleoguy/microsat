#load in packages needed
library(phytools)
library(geiger)
library(phylolm)

#read in nmicrosatellite data
str <- read.csv("../data/traits/micro.vs.chrom.csv")

#read in chromosome data
chrom <- read.csv("../data/traits/data.invert.csv")

#vector of names from chromosome data
chrom.names <- paste(chrom$Genus, chrom$species)

#loop that finds species diploid chromosome number from our data frame
for(i in 1:nrow(str)){
  # if species in microcounter matches one in chromosome data 
  if(str$species[i] %in% chrom.names){
    #store the name in vector hit
    hit <- which(chrom.names == str$species[i])[1]
    #fill in diploid number for those species that have a match in the 
    #chromosome data
    str$diploid.num[i] <- chrom$Chromosome.number..female..2N[hit]
  }
}

#make a new data frame that will contain only those species with both 
#microsatellite and chromosome number data
dat.intersect <- str[complete.cases(str),]

#write the CSV with microsatellite, genome, and diploid chromosome number data
write.csv(dat.intersect, "dat.intersect.csv")

#have to manually put in _ between species names and remove first column
#for following code to work

#then reread in the data
dat.intersect <- read.csv("../figures/dat.intersect.csv", 
                          as.is = T,
                          row.names = 1)

#read in trees
trees <- read.nexus("../data/trees/post.nex")

#make a vector to store p-values
pvals.chrom <- c()

#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
for(i in 1:100){
  #stores tree number
  tree.test <- trees[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=dat.intersect)
  #stores current trees data
  tree.cur <- foo[[1]]
  #creates data frame of the data for each tree
  dat <- as.data.frame(foo[[2]])
  #stores p-value on phylolm analysis
  pvals.chrom[i] <- summary(phylolm((all/gsz) ~ diploid.num, 
                                   data = dat.intersect, 
                                   phy = tree.cur, 
                                   model = "BM", 
                                   boot = 100))$coefficients[2,6]
}

#makes a histogram containing the p-values from the loop
hist(pvals.chrom,
     main = "Chromosome Number and Microsatellites P-Values",
     xlab = "P-Values",
     ylab = "Frequency of P-Values")

#store the bp/Mbp microsatellite content
bpMbp <- dat.intersect$all/(dat.intersect$gsz/1000000)

#plot the microsatellite content in bp/Mbp and the diploid chromosome number
plot(bpMbp~dat.intersect$diploid.num,
     xlab = "Diploid Chromosome Number",
     ylab = "Microsatellite Content (bp/Mbp)",
     pch = 16, 
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

# export pdf at 4.3" x 4.3"