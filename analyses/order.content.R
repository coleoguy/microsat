#load in necessary packages
library(geiger)

#read in insect phylogeny
trees <- read.nexus("../data/post.nex")

#read in the microsatellite data
dat.mic <- read.csv("../results/micRocounter_results_TII.csv", 
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

# make named vector for orders
order <- as.factor(dat.mic$order)
names(order) <- row.names(dat.mic)

#create results data frame and indicate proper column names
results.Mbp <- matrix(NA, 100, 2)
colnames(results.Mbp) <- c("wophylo","wphylo")

#run phyloANOVA for bpMbp and orders
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
  results.Mbp[i, 1] <- aov.sum$`Pr(>F)`[1]
  results.Mbp[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

#save p-value data frame into csv
write.csv(results.Mbp, "../results/Mbp.csv")

#make named vector for twomers
twomers <- dat.mic$twomers
names(twomers) <- row.names(dat.mic)

#create results data frame and indicate proper column names
results.twomers <- matrix(NA, 100, 2)
colnames(results.twomers) <- c("wophylo","wphylo")

#run phyloANOVA for twomers and orders
twomers.2 <- twomers
ord2 <- order
for(i in 1:100){
  for(j in 1:length(twomers.2)){
    hit <- which(names(twomers) == trees.pruned[[i]]$tip.label[j])
    twomers.2[j] <- twomers[hit]
    ord2[j] <- order[hit]
  }
  names(twomers.2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(twomers.2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results.twomers[i, 1] <- aov.sum$`Pr(>F)`[1]
  results.twomers[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

#save p-value data frame into csv
write.csv(results.Mbp, "../results/twomers.csv")

#make named vector for threemers
threemers <- dat.mic$threemers
names(threemers) <- row.names(dat.mic)

#create results data frame and indicate proper column names
results.threemers <- matrix(NA, 100, 2)
colnames(results.threemers) <- c("wophylo","wphylo")

#run phyloANOVA for twomers and orders
threemers.2 <- threemers
ord2 <- order
for(i in 1:100){
  for(j in 1:length(threemers.2)){
    hit <- which(names(threemers) == trees.pruned[[i]]$tip.label[j])
    threemers.2[j] <- threemers[hit]
    ord2[j] <- order[hit]
  }
  names(threemers.2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(threemers.2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results.threemers[i, 1] <- aov.sum$`Pr(>F)`[1]
  results.threemers[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

#save p-value data frame into csv
write.csv(results.Mbp, "../results/threemers.csv")

#make named vector for threemers
fourmers <- dat.mic$fourmers
names(fourmers) <- row.names(dat.mic)

#create results data frame and indicate proper column names
results.fourmers <- matrix(NA, 100, 2)
colnames(results.fourmers) <- c("wophylo","wphylo")

#run phyloANOVA for twomers and orders
fourmers.2 <- fourmers
ord2 <- order
for(i in 1:100){
  for(j in 1:length(fourmers.2)){
    hit <- which(names(fourmers) == trees.pruned[[i]]$tip.label[j])
    fourmers.2[j] <- fourmers[hit]
    ord2[j] <- order[hit]
  }
  names(fourmers.2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(fourmers.2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results.fourmers[i, 1] <- aov.sum$`Pr(>F)`[1]
  results.fourmers[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

#save p-value data frame into csv
write.csv(results.Mbp, "../results/fourmers.csv")

#make named vector for threemers
fivemers <- dat.mic$fivemers
names(fivemers) <- row.names(dat.mic)

#create results data frame and indicate proper column names
results.fivemers <- matrix(NA, 100, 2)
colnames(results.fivemers) <- c("wophylo","wphylo")

#run phyloANOVA for twomers and orders
fivemers.2 <- fivemers
ord2 <- order
for(i in 1:100){
  for(j in 1:length(fivemers.2)){
    hit <- which(names(fivemers) == trees.pruned[[i]]$tip.label[j])
    fivemers.2[j] <- fivemers[hit]
    ord2[j] <- order[hit]
  }
  names(fivemers.2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(fivemers.2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results.fivemers[i, 1] <- aov.sum$`Pr(>F)`[1]
  results.fivemers[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}
#save p-value data frame into csv
write.csv(results.Mbp, "../results/fivemers.csv")

#make named vector for threemers
sixmers <- dat.mic$sixmers
names(sixmers) <- row.names(dat.mic)

#create results data frame and indicate proper column names
results.sixmers <- matrix(NA, 100, 2)
colnames(results.sixmers) <- c("wophylo","wphylo")

#run phyloANOVA for twomers and orders
sixmers.2 <- sixmers
ord2 <- order
for(i in 1:100){
  for(j in 1:length(sixmers.2)){
    hit <- which(names(sixmers) == trees.pruned[[i]]$tip.label[j])
    sixmers.2[j] <- fourmers[hit]
    ord2[j] <- order[hit]
  }
  names(sixmers.2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(sixmers.2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results.sixmers[i, 1] <- aov.sum$`Pr(>F)`[1]
  results.sixmers[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

#save p-value data frame into csv
write.csv(results.Mbp, "../results/sixmers.csv")

#make named vector for threemers
all <- dat.mic$all
names(all) <- row.names(dat.mic)

#create results data frame and indicate proper column names
results.all <- matrix(NA, 100, 2)
colnames(results.all) <- c("wophylo","wphylo")

#run phyloANOVA for twomers and orders
all.2 <- all
ord2 <- order
for(i in 1:100){
  for(j in 1:length(all.2)){
    hit <- which(names(all) == trees.pruned[[i]]$tip.label[j])
    all.2[j] <- all[hit]
    ord2[j] <- order[hit]
  }
  names(all.2) <- names(ord2) <- trees.pruned[[i]]$tip.label
  fit <- aov.phylo(all.2~ord2,
                   phy = trees.pruned[[i]],
                   nsim = 100)
  aov.sum <- attributes(fit)$summary
  results.all[i, 1] <- aov.sum$`Pr(>F)`[1]
  results.all[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

#save p-value data frame into csv
write.csv(results.Mbp, "../results/all.csv")
