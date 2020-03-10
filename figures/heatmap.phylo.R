#load in packages that are needed
library(phytools)

#load in the trees
trees <- read.nexus("data/tree/tree.nex")
#select a single tree to use
tree <- trees[[sample(1:100, 1)]]
#load in the data
dat.microsat <- read.csv("results/micRocounter_results_TII.csv")

# drops the tip for B.terrestis
pruned.tree <- drop.tip(phy=tree, tip="B.terrestris")
pruned.tree <- drop.tip(phy=pruned.tree, tip="Plutella_xylostella")

#subset


phylo.heatmap(pruned.tree, fsize = c(0.10, 0.5, 1), log(dat.microsat[5:10]), hcl.colors(n= 500, palette= "viridis"))

              