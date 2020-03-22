#set figures as working directory
#load in libraries
library(phytools)
library(geiger)
library(phylolm)

#read in csv with rates of evolution
rates <- read.csv("../analyses/tip.rates/tip.rates.csv",
                  row.names = 1)

#store the average rate in a named vector by species name
rates.species <- rates$Average
names(rates.species) <- row.names(rates)

#load in chromosomes number data
dat.intersect <- read.csv("../data/traits/dat.intersect.csv",
                          as.is = T, 
                          row.names = 1)

#make an empty column in the data frame with diploid chromosome number
dat.rates.chrom <- cbind(dat.intersect, rates.evol = "", stringsAsFactors = F)
#loop that finds species diploid chromosome number from our data frame
for(i in 1:nrow(dat.rates.chrom)){
  # if species in microcounter matches one in chromosome data 
  if(row.names(dat.rates.chrom)[i] %in% names(rates.species)){
    #store the name in vector hit
    hit <- which(names(rates.species) == row.names(dat.rates.chrom)[i])
    #fill in rates for those species that have a match in the 
    #chromosome data
    dat.rates.chrom$rates.evol[[i]] <- rates.species[[hit]]
  }
}

#clean environment
rm(list = c("dat.intersect", "rates", "i", "rates.species"))

#read in trees
trees <- read.nexus("../data/trees/post.nex")

#make a vector to store p-values 
pvals.rates <- c()

#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
for(i in 1:100){
  #stores tree number
  tree.test <- trees[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=dat.rates.chrom)
  #stores current trees data
  tree.cur <- foo[[1]]
  #creates data frame of the data for each tree
  dat <- as.data.frame(foo[[2]])
  #stores p-value on phylolm analysis
  pvals.rates[i] <- summary(phylolm(as.numeric(rates.evol) ~ diploid.num, 
                                    data = dat.rates.chrom, 
                                    phy = tree.cur, 
                                    model = "BM", 
                                    boot = 100))$coefficients[2,6]
}

#makes a histogram containing the p-values from the loop
hist(pvals.rates,
     main = "Chromosome Number and Rate P-Values",
     xlab = "P-Values",
     ylab = "Frequency of P-Values",
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

#export pdf as 6" x 6"

#plot the microsatellite content in bp/Mbp and the diploid chromosome number
plot(dat.rates.chrom$rates.evol~dat.rates.chrom$diploid.num,
     xlab = "Diploid Chromosome Number",
     ylab = "Microsatellite Evolution Rates",
     pch = 16, 
     col = rgb(250, 159, 181, 100,
               maxColorValue = 255))

# export pdf at 4.3" x 4.3"