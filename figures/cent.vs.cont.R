#set working directory as figures
#load in libraries
library(geiger)

# read in microsatellite and centromere data
microsat.cent <- read.csv("../results/ssr.inference/micRocounter_results_TII_typecentromere.csv", 
                          row.names = 4)

#read in centromere type data
holo.or.mono <- microsat.cent$holo.or.mono
names(holo.or.mono) <- row.names(microsat.cent)

#import trees
trees <- read.nexus("../data/trees/post.nex")

# drops the tip 
pruned.tree <- c()
for(i in 1:100){
  pruned.tree[[i]] <- drop.tip(phy=trees[[i]], tip=c("B.terrestris",
                                                     "Plutella_xylostella",
                                                     "Timema_cristinae"))
}

# make named vector for bpMbp coontent
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
           col = rgb(250, 159, 181, 100,
                     maxColorValue = 255))


# export pdf at 4.3" x 4.3"

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

#create data frame with all pvalues
pvals.aovphylo <- data.frame(pval.2mer, pval.3mer, pval.4mer, pval.5mer, pval.6mer, pval.all, pval.bpMbp)
write.csv(pvals.aovphylo, "../results/cent.type/pvals.aovphylo.csv")
