#load in necessary packages
library(phytools)

#read in insect phylogeny
trees <- read.nexus("../../data/trees/post.nex")

#read in the microsatellite data
dat.mic <- read.csv("../../results/ssr.inference/micRocounter_results_TII.csv", 
                    as.is = T, row.names = 4)


# match up tree and data

trees.pruned <- list()
#loop through to drop any unmatching data or tree tips
for(i in 1:100){
  trees.pruned[[i]] <- treedata(phy = trees[[i]], data=dat.mic)[[1]]
}
rm(trees, i)
