# READ IN DATA NEEDED FOR ANALYSES

#set figures as wd
#load phytools
library(phytools)

# read in microsatellite and centromere data
microsat.cent <- read.csv("data/micRocounter_results_TII_typecentromere.csv", 
                          row.names = 4)

# read in sampled trees
trees <- read.nexus("data/tree/tree.nex")

# select a tree at random for making a continuous trait map
tree <- trees[[sample(1:100, 1)]]
pruned.tree <- drop.tip(phy=tree, tip="B.terrestris")
pruned.tree <- drop.tip(phy=pruned.tree, tip="Plutella_xylostella")

# RUN PHYLOANOVA ANALYSES

# load in geiger library for aovphylo
library(geiger)

# run aovphylo with phylogenetic correction
# make named vector for bpMbp coontent
bp.Mbp <- microsat.cent$bp.Mbp
names(bp.Mbp) <- row.names(microsat.cent)
# make named vector for type of centromere
holo.or.mono <- microsat.cent$holo.or.mono
names(holo.or.mono) <- row.names(microsat.cent)
#run phyloANOVA for bpMbp and centromere type
aovphylo.bpMbp <- aov.phylo(bp.Mbp ~ holo.or.mono,
                     phy = tree,
                     nsim = 100)

# make named vector for all microsat content
bp.all <- microsat.cent$all
names(bp.all) <- row.names(microsat.cent)
# run phyloANOVA for all microsat content and centromere type
aovphylo.bpall <- aov.phylo(bp.all ~ holo.or.mono,
                            phy = tree,
                            nsim = 100)

# plot for presentation
boxplot(log(bp.all) ~ holo.or.mono,
        data = microsat.cent,
        outpch = NA,
        xlab = "Type of Centromere",
        ylab = "log Microsatellite Content (bp)")
stripchart(log(microsat.cent$all) ~ microsat.cent$holo.or.mono, 
           vertical = TRUE, 
           data = microsat.cent, 
           method = "jitter", 
           add = TRUE, 
           pch = 20, 
           col = 'blue',
           bg = "bisque")

#save pdf of image 5" x 5"

#make named vector for 2mer content
bp.2 <- microsat.cent$twomers
names(bp.2) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp2 <- aov.phylo(bp.2 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

#make named vector for 3mer content
bp.3 <- microsat.cent$threemers
names(bp.3) <- row.names(microsat.cent)
# run phyloANOVA for 3mers and centromere type
aovphylo.bp3 <- aov.phylo(bp.3 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

# make named vector for 4mer content
bp.4 <- microsat.cent$fourmers
names(bp.4) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp4 <- aov.phylo(bp.4 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

# make named vector for 5mer content
bp.5 <- microsat.cent$fivemers
names(bp.5) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp5 <- aov.phylo(bp.5 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)

# make named vector for 6mer content
bp.6 <- microsat.cent$sixmers
names(bp.6) <- row.names(microsat.cent)
# run phyloANOVA for 2mers and centromere type
aovphylo.bp2 <- aov.phylo(bp.6 ~ holo.or.mono,
                          phy = tree,
                          nsim = 100)
