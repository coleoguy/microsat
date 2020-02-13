#read in necessarylibraries
library(phytools)
library(geiger)
library(phylolm)

#read in the file containing the microsat and diploid number data
chroms <- read.csv("data/micRocounter/TAGresults.TIIremoved.csv", as.is = T, 
                   row.names = 1)

#read in trees
tree <- read.nexus("data/tree/tree.nex")

# make a vector ready for TAG data values
p.vals.TAG <- c()

#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
for(i in 1:100){
  #stores tree number
  tree.test <- tree[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=chroms)
  #stores current trees data
  tree.cur <- foo[[1]]
  #creates data frame of the data for each tree
  dat <- as.data.frame(foo[[2]])
  #stores p-value on phylolm analysis
  p.vals.TAG[i] <- summary(phylolm(all ~ Total.Bases..TAG, 
                               data = dat, 
                               phy = tree.cur, 
                               model = "BM", 
                               boot = 100))$coefficients[2,6]
}

#makes a histogram containing the p-values from the loop
hist(p.vals.TAG,
     main = "TAG P-Values",
     xlab = "P-Values",
     ylab = "Frequency of P-Values")

#save image as PDF 4" x 6"

# make a vector ready for TTAGG data values
p.vals.TTAGG <- c()

#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
for(i in 1:100){
  #stores tree number
  tree.test <- tree[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=chroms)
  #stores current trees data
  tree.cur <- foo[[1]]
  #creates data frame of the data for each tree
  dat <- as.data.frame(foo[[2]])
  #stores p-value on phylolm analysis
  p.vals.TTAGG[i] <- summary(phylolm(all ~ Total.Bases.TTAGG, 
                                   data = dat, 
                                   phy = tree.cur, 
                                   model = "BM", 
                                   boot = 100))$coefficients[2,6]
}
#makes a histogram containing the p-values from the loop
hist(p.vals.TTAGG,
     main = "TTAGG P-Values",
     xlab = "P-Values",
     ylab = "Frequency of P-Values")

#save image as PDF 4" x 6"

# make a vector ready for TTAGG data values
p.vals.TCAGG <- c()
#loops through the 100 posterior distribution trees and determines the data
#necessary for the p-value
for(i in 1:100){
  #stores tree number
  tree.test <- tree[[i]]
  #matches species within the dataset and the tree
  foo <- treedata(phy = tree.test, data=chroms)
  #stores current trees data
  tree.cur <- foo[[1]]
  #creates data frame of the data for each tree
  dat <- as.data.frame(foo[[2]])
  #stores p-value on phylolm analysis
  p.vals.TCAGG[i] <- summary(phylolm(all ~ Total.Loci.TCAGG.1, 
                                     data = dat, 
                                     phy = tree.cur, 
                                     model = "BM", 
                                     boot = 100))$coefficients[2,6]
}
#makes a histogram containing the p-values from the loop
hist(p.vals.TCAGG,
     main = "TCAGG P-Values",
     xlab = "P-Values",
     ylab = "Frequency of P-Values")

#save image as PDF 4" x 6"