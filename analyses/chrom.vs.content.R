#set figures as working directory
#load in packages needed
library(phytools)
library(geiger)
library(phylolm)

#read in nmicrosatellite data
str <- read.csv("../data/micro.vs.chrom.csv")

#read in chromosome data
chrom <- read.csv("../data/data.invert.csv")

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

#read in trees
trees <- read.nexus("../data/post.nex")

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