#load in packages that are needed
library(phytools)

#load in the trees
trees <- read.nexus("../data/trees/post.nex")
#select a single tree to use
tree <- trees[[sample(1:100, 1)]]
rm(trees)
#load in the data
dat.microsat <- read.csv("../results/ssr.inference/micRocounter_results_TII.csv",
                         row.names = 4)

# drops the tip for B.terrestis
pruned.tree <- drop.tip(phy=tree, tip=c("B.terrestris",
                                        "Plutella_xylostella"))



phylo.heatmap(tree = pruned.tree, 
              fsize = c(.0001, .0001, .31), standardize=T,
              X = log(dat.microsat[,4:9]*1000), 
              labels = F, pts = F,
              colors = hcl.colors(n = 500, palette = "viridis"))
exp(-3.38)
              