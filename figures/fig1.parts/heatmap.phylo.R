#set working directory as figures/heatmap
#load in packages that are needed
library(phytools)

#load in the trees
trees <- read.nexus("../../data/trees/post.nex")
#select a single tree to use
tree <- trees[[sample(1:100, 1)]]
rm(trees)
#load in the data
dat.microsat <- read.csv("../../results/ssr.inference/micRocounter_results_TII.csv",
                         row.names = 4)

# drops the tip for B.terrestis
pruned.tree <- drop.tip(phy=tree, tip=c("B.terrestris",
                                        "Plutella_xylostella",
                                        "Timema_cristinae"))
#this loop will make sure the tip labels and the microsat data are in the
#same order
foo <- dat.microsat
sp <- c()
for(i in 1:nrow(dat.microsat)){
  hit <- which(row.names(dat.microsat)==pruned.tree$tip.label[i])
  foo[i, ] <- dat.microsat[hit, ]
  sp[i] <- row.names(dat.microsat)[hit]
}
row.names(foo) <- sp

#make the heatmap with phylogeny
phylo.heatmap(tree = pruned.tree, 
              fsize = c(.0001, .0001, .31), standardize=T,
              X = log(dat.microsat[,11:16]*1000), 
              labels = F, pts = F,
              colors = hcl.colors(n = 500, palette = "viridis"))
#export as pdf 7" x 7"              




