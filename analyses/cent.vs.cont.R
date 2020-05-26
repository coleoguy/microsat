#set working directory as figures
#load in libraries
library(geiger)

# read in microsatellite and centromere data
microsat.cent <- read.csv("../results/micRocounter_results_TII_typecentromere.csv",
                          row.names = 4)

#read in centromere type data
holo.or.mono <- microsat.cent$holo.or.mono
names(holo.or.mono) <- row.names(microsat.cent)

#import trees
trees <- read.nexus("../data/post.nex")

# drops the tip
pruned.tree <- c()
for(i in 1:100){
  pruned.tree[[i]] <- drop.tip(phy=trees[[i]], tip=c("B.terrestris",
                                                     "Plutella_xylostella",
                                                     "Timema_cristinae"))
}

# make named vector for bpMbp content
bp.Mbp <- microsat.cent$bp.Mbp
names(bp.Mbp) <- row.names(microsat.cent)

# run aovphylo with phylogenetic correction
aovphylo.bpMbp <- pval.bpMbp <- c()
for(i in 1:100){
  aovphylo.bpMbp[[i]] <- aov.phylo(bp.Mbp ~ holo.or.mono,
                          phy = pruned.tree[[i]],
                          nsim = 100)
  pval.bpMbp[[i]] <- print(attributes(aovphylo.bpMbp[[i]])$summary)[1,6]
}

# make named vector for all microsat content
bp.all <- microsat.cent$all
names(bp.all) <- row.names(microsat.cent)

# run aovphylo with phylogenetic correction
aovphylo.bpall <- pval.all <- c()
for(i in 1:100){
  aovphylo.bpall[[i]] <- aov.phylo(bp.all ~ holo.or.mono,
                                   phy = pruned.tree[[i]],
                                   nsim = 100)
  pval.all[[i]] <- print(attributes(aovphylo.bpall[[i]])$summary)[1,6]
}

#make named vector for 2mer content
bp.2 <- microsat.cent$twomers
names(bp.2) <- row.names(microsat.cent)

# run phyloANOVA for 2mers and centromere type
aovphylo.bp2 <- pval.2mer <- c()
for(i in 1:100){
  aovphylo.bp2[[i]] <- aov.phylo(bp.2 ~ holo.or.mono,
                                 phy = pruned.tree[[i]],
                                 nsim = 100)
  pval.2mer[[i]] <- print(attributes(aovphylo.bp2[[i]])$summary)[1,6]
}

#make named vector for 3mer content
bp.3 <- microsat.cent$threemers
names(bp.3) <- row.names(microsat.cent)

# run phyloANOVA for 3mers and centromere type
aovphylo.bp3 <- pval.3mer <- c()
for(i in 1:100){
  aovphylo.bp3[[i]] <- aov.phylo(bp.3 ~ holo.or.mono,
                            phy = pruned.tree[[i]],
                            nsim = 100)
  pval.3mer[[i]] <- print(attributes(aovphylo.bp3[[i]])$summary)[1,6]
}

# make named vector for 4mer content
bp.4 <- microsat.cent$fourmers
names(bp.4) <- row.names(microsat.cent)

# run phyloANOVA for 2mers and centromere type
aovphylo.bp4 <- pval.4mer <- c()
for(i in 1:100){
  aovphylo.bp4[[i]] <- aov.phylo(bp.4 ~ holo.or.mono,
                                 phy = pruned.tree[[i]],
                                 nsim = 100)
  pval.4mer[[i]] <- print(attributes(aovphylo.bp4[[i]])$summary)[1,6]
}

# make named vector for 5mer content
bp.5 <- microsat.cent$fivemers
names(bp.5) <- row.names(microsat.cent)

# run phyloANOVA for 2mers and centromere type
aovphylo.bp5 <- pval.5mer <- c()
for(i in 1:100){
aovphylo.bp5[[i]] <- aov.phylo(bp.5 ~ holo.or.mono,
                          phy = pruned.tree[[i]],
                          nsim = 100)
pval.5mer[[i]] <- print(attributes(aovphylo.bp5[[i]])$summary)[1,6]
}

# make named vector for 6mer content
bp.6 <- microsat.cent$sixmers
names(bp.6) <- row.names(microsat.cent)

# run phyloANOVA for 2mers and centromere type
aovphylo.bp6 <- pval.6mer <- c()
for(i in 1:100){
aovphylo.bp6[[i]] <- aov.phylo(bp.6 ~ holo.or.mono,
                          phy = pruned.tree[[i]],
                          nsim = 100)
pval.6mer[[i]] <- print(attributes(aovphylo.bp6[[i]])$summary)[1,6]
}

# get non-phylo pvals
attributes(aovphylo.bp2[[1]]) #0.168
attributes(aovphylo.bp3[[1]]) #0.175
attributes(aovphylo.bp4[[1]]) #0.364
attributes(aovphylo.bp5[[1]]) #0.717
attributes(aovphylo.bp6[[1]]) #0.807
attributes(aovphylo.bpall[[1]]) #0.185
attributes(aovphylo.bpMbp[[1]]) #0.109

#create data frame with all pvalues
pvals.aovphylo <- data.frame(pval.2mer, pval.3mer, pval.4mer, pval.5mer, pval.6mer, pval.all, pval.bpMbp)
write.csv(pvals.aovphylo, "../results/cent.type/pvals.aovphylo.csv")
