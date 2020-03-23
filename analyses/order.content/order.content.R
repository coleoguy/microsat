#load in necessary packages
library(geiger)

#read in insect phylogeny
trees <- read.nexus("../../data/trees/post.nex")

#read in the microsatellite data
dat.mic <- read.csv("../../results/ssr.inference/micRocounter_results_TII.csv", 
                    as.is = T, row.names = 4)

#loop through to drop any unmatching data or tree tips
trees.pruned <- c()
for(i in 1:100){
  trees.pruned[[i]] <- treedata(phy = trees[[i]], data=dat.mic)[[1]]
}

# run aovphylo with phylogenetic correction
# make named vector for bpMbp coontent
bp.Mbp <- dat.mic$bp.Mbp
names(bp.Mbp) <- row.names(dat.mic)
# make named vector for type of centromere
order <- as.factor(dat.mic$order)
names(order) <- row.names(dat.mic)
#run phyloANOVA for bpMbp and centromere type
results <- matrix(NA, 100, 2)
colnames(results) <- c("wophylo","wphylo")
for(i in 1:100){
  fit <- aov.phylo(bp.Mbp ~ order,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results[i, 1] <- aov.sum$`Pr(>F)`[1]
  results[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}
write.csv(results,file="../results/cent.vs.cont.csv",row.names = F)

