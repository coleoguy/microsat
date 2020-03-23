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


bp2 <- bp.Mbp
ord2 <- order
for(i in 1:100){
  for(j in 1:length(bp2)){
    hit <- which(names(bp.Mbp) == trees.pruned[[i]]$tip.label[j])
    bp2[j] <- bp.Mbp[hit]
    ord2[j] <- order[hit]
  }
  names(bp2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(bp2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results[i, 1] <- aov.sum$`Pr(>F)`[1]
  results[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

twomers <- dat.mic$twomers
names(twomers) <- row.names(dat.mic)
twomers.2 <- twomers
ord2 <- order
for(i in 1:100){
  for(j in 1:length(twomers.2)){
    hit <- which(names(twomers) == trees.pruned[[i]]$tip.label[j])
    twomers.2[j] <- twomers[hit]
    ord2[j] <- order[hit]
  }
  names(bp2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(twomers.2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results[i, 1] <- aov.sum$`Pr(>F)`[1]
  results[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}



write.csv(results,file="../results/cent.vs.cont.csv",row.names = F)

