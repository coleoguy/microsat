#set working directory as figures/heatmap
#load in packages that are needed
library(phytools)

#load in the trees
trees <- read.nexus("../data/post.nex")
#select a single tree to use
tree <- trees[[sample(1:100, 1)]]
rm(trees)

#load in tip rate data
tips <- read.csv("../results/tip.rates.csv", row.names = 1)
#subset out columns that we need
tips <- tips[101]

#this loop will make sure the tip labels and the microsat data are in the
#same order
foo2 <- tips
sp2 <- c()
for(i in 1:nrow(tips)){
  hit2 <- which(row.names(tips)==pruned.tree$tip.label[i])
  foo2[i, ] <- tips[hit2, ]
  sp2[i] <- row.names(tips)[hit2]
}
row.names(foo2) <- sp2

tips <- foo2$Average
names(tips) <- row.names(foo2)

#plot tree with bars
plotTree.wBars(tree = pruned.tree,
               x = abs(tips))

#export as pdf 7" x 7" 