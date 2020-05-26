# load in geiger library for aovphylo
library(geiger)

# read in microsatellite and centromere data
microsat.cent <- read.csv("../results/micRocounter_results_TII_typecentromere.csv", 
                          row.names = 4)

# read in sampled trees
trees <- read.nexus("../data/post.nex")

# RUN PHYLOANOVA ANALYSES

# run aovphylo with phylogenetic correction
# make named vector for bpMbp coontent
bp.Mbp <- microsat.cent$bp.Mbp
names(bp.Mbp) <- row.names(microsat.cent)

# make named vector for type of centromere
holo.or.mono <- microsat.cent$holo.or.mono
names(holo.or.mono) <- row.names(microsat.cent)

#run phyloANOVA for bpMbp and centromere type
results <- matrix(NA, 100, 2)
colnames(results) <- c("wophylo","wphylo")
for(i in 1:100){
  fit <- aov.phylo(bp.Mbp ~ holo.or.mono,
                              phy = trees[[i]],
                              nsim = 100)
  aov.sum <- attributes(fit)$summary
  results[i, 1] <- aov.sum$`Pr(>F)`[1]
  results[i, 2] <- aov.sum$`Pr(>F) given phy`[1]
}

#write results to a csv
write.csv(results,file="../results/cent.vs.cont.csv",row.names = F)



