#load phytoools
library(phytools)

#read in the first 50 trees from BEAST run 1
# RENAME SAMPLED TREES
tree1 <- read.nexus("pruned .nex")

#read in the second 50 trees from BEAST run 2
tree2 <- read.nexus("pruned(2) .nex")

#get all 100 trees into one object
trees <- c(tree1, tree2)
rm(tree1,tree2)

#prepare tree for making a continuous trait map
tree <- trees[[sample(1:100, 1)]]

#read in  microsatellite data 
dat.microsat <- read.csv("../../../results/bp.results.TII.csv", as.is=T,
                         row.names = 4)


library(geiger)

# match up tree and data
# make a vector ready for treedata
mscont <- dat.microsat[,18]
names(mscont) <- row.names(dat.microsat)
tree <- treedata(phy = tree, data=mscont)[[1]]

write.tree(tree, file = "tree.BAMM")


library(ape)
BAMM.tree <- read.tree("tree.BAMM.tre")
plot(BAMM.tree)

#tree needs to be ultrametric and binary
is.ultrametric(BAMM.tree)
is.binary(BAMM.tree)


